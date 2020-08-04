import os, sys
import subprocess
import re, glob
import time

## Assumptions/naming conventions:

## all raw data for a sample is in a specified (derived from sampleNr) sample dir
## raw data is provided as fastq.gz files somewhere in that directory
## forward and reverse read files are named identical except _R1_ is _R2_ in the reverse
## scratch partition is writable for script owner and has enough space
## assuming trimmomatic produces gzipped fastq files! (*.fastq.gz)
## Lane information is in the file name like this: ..._L001_...

## Reference genome must be indexed (bwa -index)


def extractLane(fname): return re.findall('_L\d\d\d_', os.path.basename(fname))[-1].strip()[-2]

## Global variable setting
if True:
    outdir = "/research/btc_bioinformatic/operations/scratch_WES/"
    tmpdir = "/research/btc_bioinformatic/operations/scratch_WES/tmp/"
    
    applicationDir = '/research/btc_bioinformatic/operations/'
    rawData = '/research/btc_bioinformatic/results/'
    resourceDir = '/research/btc_bioinformatic/operations/data_masdar/'
    Mills = "%s/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf" % resourceDir
    G1000_indels = "%s/1000G_phase1.indels.hg19.sites.vcf"% resourceDir
    G1000_snps = "%s/1000G_phase1.snps.high_confidence.hg19.sites.vcf"% resourceDir
    Hapmap = "%s/hapmap_3.3.hg19.sites.vcf"% resourceDir
    Omni = "%s/1000G_omni2.5.hg19.sites.vcf"% resourceDir
    dbsnp = "%s/dbsnp_138.hg19.vcf"% resourceDir  # we should change to a most recent dbsnp
    TR = "%s/truseq-exome-targeted-regions-manifest-v1-2a.bed"% resourceDir  # Target region for sequencing
    jobOutputDir = '%s/Pipeline/JobOutputs' % applicationDir

    fastqcDir = "%s/FastQC" % applicationDir
    picardDir = '%s/picard/build/libs' % applicationDir
    TrimmomaticsDir = "%s/Trimmomatic-0.36" % applicationDir
    illuminaAdapterFile = '%s/AuxData/TruSeq3-PE.fa' % applicationDir
    qualimapDir = "%s/qualimap_v2.2.1" % applicationDir
    rawDataDir = '%s/WES_NF' % rawData

    refDirLookup = {'hg19':'%s/data_masdar/ucsc.hg19.fasta' % applicationDir, 
                    'hg38':'%s/AuxData/hg38/Homo_sapiens_assembly38.fasta' % applicationDir}

pipelineOrder = ['trimming', 'fastqc', 'bwa', 'SortSam', 'MergeBam', 'MarkDuplicates',
                 'qualimap', 'BuildBamIndex', 'BaseRecalibrator', 'ApplyBQSR', 'HaplotypeCaller']
class Pipeline:
    ## static class variables

    ## TODO check out spark version
    ## Introduce Checkpoints for automated pipeline reruns

    def __init__(self, sampleNr=0, refShort='hg19', singleLane=True, cleanup='asYouGo', sample=None, picardCmd='local'):
        self.cleanup = cleanup ## Possible str values: asYouGo, atTheEnd, none
        self.refShort = refShort
        self.ref=refDirLookup[refShort]
        self.sampleNr = sampleNr
        self.singleLane = singleLane
        self.completedJobs = []

        if singleLane: ## Second Batch, different naming conventions
            if sample is None:
                self.rawDataDir = sorted(glob.glob('%s/1*'%rawDataDir))[sampleNr-1]
                self.sample = os.path.basename(self.rawDataDir)
            else:
                self.sample = sample
                self.rawDataDir = '%s/%s'%(rawDataDir,sample)
        else: ## Local (Mariam's) samples
            if not sample is None:
                self.rawDataDir = '/research/btc_bioinformatic/results/WES/%s' % sample
                self.sample = sample

        print("Sample: %s"%self.sample)
        sys.stdout.flush()
        self.outdir = '%s/%s' % (outdir, self.sample)
        self.dirs = {'outdir' : self.outdir,
                     'initfastqOutdir' : '%s/FastQCinit'%self.outdir,
                     'fastqOutdir' : '%s/FastQC'%self.outdir,
                     'qualimapOutdir' : '%s/Qualimap'%self.outdir,
                     'trimOutdir' : '%s/Trim'%self.outdir,
                     'alignOutdir' : '%s/Alignment'%self.outdir,
                     'gatkOutdir' : '%s/GATK'%self.outdir}

        for location in self.dirs.values():
            if not os.path.exists(location): os.mkdir(location)

        self.picardCmdLookup = {'local':picardDir + '/picard.jar',
                                'environment': '$PICARD'}
        self.picardCmd = self.picardCmdLookup.get(picardCmd, picardCmd)
        self.zipped = ''

    def findFiles(self, inputDir, pattern="*.fastq.gz"):
        fastqFiles0 = subprocess.check_output(["find",  inputDir, "-name", pattern ])
        if len(fastqFiles0) == 0:
            print ("Warning, %s has no files with pattern '%s'!" % (inputDir, pattern))
            sys.stdout.flush()
        return sorted([ff for ff in fastqFiles0.decode().split('\n') if ff])
    def determineSampleID(self, dir):
        ## scraping the sample ID out of the subdirectory - relies heavily on naming convention!
        fn = self.findFiles(dir)[0] ## not efficient!
        return fn[len(dir):].split('/')[2].split('-')[0]
    def fastQC(self, inputDir, outputDir, blocking=False):
        ## TODO: refactor against new Job class
        fastqFiles = self.findFiles(inputDir)
        print ("Dealing with %s fastqFiles" % len(fastqFiles))
        if len(fastqFiles) == 0:
            print ("Warning, %s has no fastq.gz files!" % inputDir)
            return
        threads = int(24/len(fastqFiles))
        processes = []
        for ff in fastqFiles[:1]:   ## CHANGE!!!!
            cmd = ['%s/%s' % (fastqcDir, 'fastqc'),'--threads=%d' % threads, '--outdir=%s'%outputDir, ff]
            print(" ".join(cmd))
            processes.append(subprocess.Popen(cmd))
        print("Finished fastqc, results in %s"% outputDir)
        if blocking: ## This can run independently
            for p in processes:
                p.wait()
            print ("[PIPE:] Finished FastQC")
        sys.stdout.flush()

    ## *fastq.gz -> _*[un]paired.fastqz
    def trim(self, pattern="*.fastq.gz"):
        fastqFiles = self.findFiles(self.rawDataDir, pattern=pattern)
        jobs = []
        for fastqFile, reverse in zip(fastqFiles[::2], fastqFiles[1::2]): ##  assume every 2. file is resp. reverse
            outBasename    = os.path.basename(fastqFile)[:-9]
            outRevBasename = os.path.basename(reverse)[:-9]
            trimmomaticCmd = ['java', '-Xmx32G', '-jar', '%s/trimmomatic-0.36.jar' % TrimmomaticsDir,
                              'PE', '-threads', '30', '-phred33',
                              fastqFile, reverse,
                              '%s/%s.paired.fastq.gz' % (self.dirs['trimOutdir'], outBasename),
                              '%s/%s.unpaired.fastq.gz' % (self.dirs['trimOutdir'], outBasename),
                              '%s/%s.paired.fastq.gz' % (self.dirs['trimOutdir'], outRevBasename),
                              '%s/%s.unpaired.fastq.gz' % (self.dirs['trimOutdir'], outRevBasename),
                              'ILLUMINACLIP:%s:2:30:10' % illuminaAdapterFile,
                              'CROP:74','SLIDINGWINDOW:4:20','MINLEN:36']
            job = Job(trimmomaticCmd, run=False)
            job.startProcess()
            jobs.append(job)
        ## Wait til processes are finished
        for job in jobs: job.wait()
        self.completedJobs.append(jobs)

    ##paired.fastqz -> paired_aligned.sam
    def align(self, pattern="*.paired.fastq.gz", cleanup=True):
        ## TODO: mapping unpaired?  if they make up for more than say 8%
        fastqFiles = self.findFiles(self.dirs['trimOutdir'], pattern=pattern)
        jobs = []
        for fastqFile, reverse in zip(fastqFiles[::2], fastqFiles[1::2]):
            if self.singleLane: lane='001'
            else: lane = extractLane(fastqFile)
            samfilename = '%s/%s_%s_FR_%s.paired_aligned.sam' % (self.dirs['alignOutdir'], self.sample, lane, self.refShort)
            alignCmd = 'bwa mem -M -t 46 -R ' +\
                        r'"@RG\tID:ID.%s\tSM:%s\tPL:Illumina\tLB:KU_WES\tPU:PU.%s" ' % (lane, self.sample, lane) +\
                        '%s %s %s'  % (self.ref, fastqFile, reverse)
            
            job = Job([alignCmd], run=False)
            if cleanup:
                job.mark4cleanup(fastqFile, condition=samfilename)
                job.mark4cleanup(reverse, condition=samfilename) 
            job.startProcess(stdout=samfilename)
            jobs.append(job)

        for j in jobs: j.wait()
        self.completedJobs.append(jobs)

    ## paired_aligned.sam ->  Sorted.bam
    def sortSam(self, pattern="*.sam", blocking=False, jvmSpace=32,
                jvmSpaceDivide=True):
        ## Requires picard module (that sets $PICARD)
        jobs = []
        samFiles = self.findFiles(self.dirs['alignOutdir'], pattern=pattern)
        if jvmSpaceDivide: jvmSpace = int(jvmSpace / (len(samFiles)))
        for samFile in samFiles:
            sortedBamFile = samFile[:-4] + '_Sorted.bam'  ## prob. nicer with os.splitext
            sortCmd = ['java', '-Xmx%sG' % jvmSpace, '-jar', self.picardCmd, 'SortSam', 'CREATE_INDEX=true',
                       'I=%s' % samFile, 'O=%s' % sortedBamFile, 'SORT_ORDER=coordinate',
                       'VALIDATION_STRINGENCY=STRICT', 'TMP_DIR=%s' % tmpdir]
            job = Job(sortCmd, run=False)
            job.mark4cleanup(samFile, condition=sortedBamFile)
            job.startProcess()
            if blocking: job.wait()
            else: jobs.append(job)
        for j in jobs: j.wait() ## only loops if there are non-blocking jobs
        self.completedJobs.append(jobs)

    ## Sorted.bam -> Sorted_merged.bam
    def mergeBam(self, jvmSpace=32):
        bamFiles = self.findFiles(self.dirs['alignOutdir'], pattern="*_Sorted.bam")
        if len(bamFiles) == 1: ## if single Lane, simply rename/link so to keep name convention for next steps
            mergedBamFile = bamFiles[0][:-4] + '_Merged.bam' ## prob. nicer with os.splitext
            subprocess.Popen('mv %s %s' % (bamFiles[0], mergedBamFile), shell=True).wait() ## ln -s seems to be slow
            return []
        mergeCmd = ['java', '-Xmx32G', '-jar', self.picardCmd, 'MergeSamFiles', 'USE_THREADING=true']
        for bamFile in bamFiles:
            mergeCmd.append('I=%s'%bamFile)
        mergedBamFile = bamFiles[0][:-4] + '_Merged.bam' ## prob. nicer with os.splitext
        mergeCmd.append('O=%s'%mergedBamFile)
        tmpfiles = '%s/*_Sorted.bam' % self.dirs['alignOutdir']
        self.completedJobs.append(Job(mergeCmd, mark4del=(tmpfiles, mergedBamFile)))

    ## Sorted_merged.bam -> Sorted_merged_dedup.bam ## this bam should be kept!
    def markDuplicates(self, jvmSpace=32):
        bamFile = self.findFiles(self.dirs['alignOutdir'], pattern="*_Sorted_Merged.bam")[0]
        
        nrBamFile = bamFile[:-4] + '_dedup.bam' ## prob. nicer with os.splitext
        metricsFile = bamFile[:-4] + '_metrics.txt' 
        dupCmd = ['java', '-Xmx%sG'%jvmSpace, '-jar', self.picardCmd, 'MarkDuplicates',
                  'I=%s'%bamFile, 'O=%s'%nrBamFile, 'REMOVE_DUPLICATES=true', 'M=%s'%metricsFile]
        self.completedJobs.append(Job(dupCmd, mark4del=(bamFile, nrBamFile)))

    def qualimap(self, blocking=False):
        bamFile = self.findFiles(self.dirs['alignOutdir'], pattern="*_dedup.bam")[0]
        qmCmd = ['%s/qualimap' % qualimapDir, 'bamqc', '-bam', bamFile, '-outdir', self.dirs['qualimapOutdir'],'-gff',TR, '--java-mem-size=32G']
        process = subprocess.Popen(" ".join(qmCmd), shell=True)
        print ("[PIPE:] %s" " ".join(qmCmd))
        if blocking:
            process.wait()
            print ("[PIPE:] Finished Qualimap")
            
    def buildBamIndex(self, jvmSpace=32): 
        bamFile = self.findFiles(self.dirs['alignOutdir'], pattern="*_dedup.bam")[0]
        bbiCmd = ['java', '-Xmx%sG'%jvmSpace, '-jar', self.picardCmd, 'BuildBamIndex', 'I=%s'%bamFile]
        self.completedJobs.append(Job(bbiCmd))

    ## Sorted_merged_dedup.bam -> Sorted_merged_dedup_recal_table.grp
    def baseRecalibrator(self, jvmSpace=32):
        bamFile = self.findFiles(self.dirs['alignOutdir'], pattern="*_dedup.bam")[0]
        #output = bamFile[:-4] + '_realigner.intervals'
        output = bamFile[:-4] + '_recal_table.grp'
        recCmd = ['gatk', '--java-options', '"-Xmx%sG"'%jvmSpace, 'BaseRecalibrator',
                  '-R', self.ref, '-I', bamFile, '-L', TR,
                  '--known-sites', G1000_snps,
                  '--known-sites', Mills,
                  '-O', output]
        self.completedJobs.append(Job(recCmd, run=True))

    ## Sorted_merged_dedup.bam -> Sorted_merged_dedup_BQSR.bam (huge bam, delete when possible!)
    def applyBQSR(self, jvmSpace=32):
        bamFile = self.findFiles(self.dirs['alignOutdir'], pattern="*_dedup.bam")[0]
        table = bamFile[:-4] + '_recal_table.grp' ## assume table produced as in baseRecalibrator
        output = bamFile[:-4] + '_BQSR.bam'
        recCmd = ['gatk', '--java-options', '"-Xmx%sG"'%jvmSpace, 'ApplyBQSR',
                  '-R', self.ref,
                  '-I', bamFile,
                  '--bqsr-recal-file', table,
                  '-O', output]
        self.completedJobs.append(Job(recCmd, run=True))

    ## Sorted_merged_dedup_BQSR.bam ->  $GATKdir/<sample>_HaplotypeCaller.vcf
    def haplotypeCaller(self, jvmSpace=32, filetype='gvcf', compress=True):
        self.filetype = filetype ## remember to set the filetype if haplotypeCaller is bypassed/omitted!!
        bamFile = self.findFiles(self.dirs['alignOutdir'], pattern="*_BQSR.bam")[0]
        output = '%s/%s_HaplotypeCaller.%s' %(self.dirs['gatkOutdir'], self.sample, filetype)
        if compress:
            output += '.gz'
            self.zipped = '.gz'

        recCmd = ['gatk', '--java-options', '"-Xmx%sG"'%jvmSpace, 'HaplotypeCaller',
                  '-R', self.ref,
                  '-I', bamFile, '-O', output,
                  '--dbsnp', dbsnp, '--genotyping-mode', 'DISCOVERY']
        if filetype=='gvcf':
            recCmd.append('-ERC GVCF') 
        self.completedJobs.append(Job(recCmd, run=True, mark4del=(bamFile, output)))

    ## Sorted_merged_dedup_BQSR.vcf ->  Sorted_merged_dedup_BQSR.{recal.vcf, tranches, rplots.R} 
    def variantRecalibrator(self, jvmSpace=32, mode='SNP', pattern="*_HaplotypeCaller"):
        vcfFile = self.findFiles(self.dirs['gatkOutdir'], pattern=pattern + '.' + self.filetype + self.zipped)[0]
        base = os.path.splitext(vcfFile)[0]
        recal_output = base + '_recal_%s.%s' % (mode, self.filetype)
        tranches_output = base + '_recal_%s.tranches' % mode
        plots_output = base + '_recal_%s.rplots.R' % mode
        G1000 = {'SNP': G1000_snps, 'INDEL': G1000_indels}[mode]
        if mode=='SNP':
            modeSpecificResources = ['--resource hapmap,known=false,training=true,truth=true,prior=15.0:%s' % Hapmap,
                                     '--resource omni,known=false,training=true,truth=true,prior=12.0:%s' % Omni]
        elif mode=='INDEL':
            modeSpecificResources = ['--resource:mills,known=false,training=true,truth=true,prior=12.0:%s' % Mills]

        recCmd = ['gatk', '--java-options', '"-Xmx%sG"' % jvmSpace, 'VariantRecalibrator',
                  '-R', self.ref,
                  '-V', vcfFile,
                  '-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum',
                  '--mode', mode,
                  '-O', recal_output, '--tranches-file', tranches_output,
                  '--rscript-file', plots_output,
                  '--resource 1000G,known=false,training=true,truth=false,prior=10.0:%s' % G1000,
                  '--resource dbsnp,known=true,training=false,truth=false,prior=2.0:%s' % dbsnp] + modeSpecificResources
        self.completedJobs.append(Job(recCmd, run=True))
        # return Code 3 when no data is found (e.g. on toy genome)

    def applyVQSR(self, jvmSpace=32, filterLevel='99.0', mode='SNP', pattern="*_HaplotypeCaller"):
        vcfFile = self.findFiles(self.dirs['gatkOutdir'], pattern=pattern + '.' + self.filetype)[0]
        base = os.path.splitext(vcfFile)[0]
        output = base + '_recal_%s.vcf' % mode
        recal_input = base + '.recal.%s' % self.filetype
        tranches_input = base + '.recal.%s' % self.filetype
        recCmd = ['gatk', '--java-options', '"-Xmx%sG"' % jvmSpace, 'ApplyVQSR',
                  '-R', self.ref,
                  '-V', vcfFile,
                  '-O', output,
                  '--ts_filter_level', filterLevel,
                  '--recal-file',  recal_input,
                  '--tranches-file', tranches_input,
                  '--mode ', mode]
        self.completedJobs.append(Job(recCmd, run=True))

    def report(self, last=-1, stopOnFail=True): #use last=0 for all
        cleanup = (self.cleanup == 'asYouGo' and last==-1) or (self.cleanup == 'atTheEnd' and last==0)
        for jobs in self.completedJobs[last:]:
            if type(jobs)!=type([]): jobs = [jobs] ## for convenience, allow singleton jobs on the list
            allOK = all([job.evaluate(cleanup=cleanup) for job in jobs])
            if not allOK:
                print ("[PIPE:] Stopping gracefully because of previous errors")
                sys.exit(-1)
            print ("[PIPE:] Finished %s successfully/%s job(s). Time %s" % (jobs[0].short, len(jobs), time.asctime()))
            print ("[PIPE:] Job finishing time(s): %s" % [job.finishingTime for job in jobs])
            print ("[PIPE:] Job duration(s) in sec: %s" % [job.duration for job in jobs])
            print ("[PIPE:] Job return code(s): %s" % [job.returnCode for job in jobs])
        sys.stdout.flush()

    def getLastCheckpoint(self):
        from restartPipe import getLastCheckpoint
        finishedJobs = getLastCheckpoint(self.sample)
        try:
            lastJob = finishedJobs[-1]

            idx = pipelineOrder.index(lastJob)
            print ("Trying to pick up after %s (job idx: %s)"% (lastJob, idx))
            return idx
        except (ValueError, IndexError) as e:
            print ("Warning, could not identify last finished job, starting from scratch!", finishedJobs)
            return -1
        sys.stdout.flush()

class Job:
    def __init__(self, job, mark4del=None, run=True, verbose=False):
        self.cmd = " ".join(job)
        if verbose:
            print ("[PIPE:] %s" % self.cmd)
        self.job = job
        self.marked4deletion = [mark4del] if mark4del else []
        if len(job) == 1: self.short=job[0].split()[0]
        elif job[0]=='java' and 'trim' in job[3]: self.short = 'trimming'
        elif job[0]=='java' and len(job)>4: self.short = job[4]
        elif job[0]=='gatk': self.short = 'GATK ' + job[3]
        else: self.short = self.cmd[:15]

        if run:
            self.startProcess()
            self.wait()
    def startProcess(self, stdout=None, shell=True):
        self.startTime = time.time()
        if stdout:
            self.cmd += ' > %s' % stdout
            with open(stdout, "w") as fh:
                self.process = subprocess.Popen(self.cmd, shell=shell, stdout=fh)
        else:
            print ("Starting process %s" % self.cmd)
            self.process = subprocess.Popen(self.cmd, shell=shell)

    def wait(self):
        self.returnCode = self.process.wait()
        self.duration = pptime (time.time() - self.startTime)
        self.finishingTime = time.asctime()
    def evaluate(self, verbose=False, cleanup=True):
        if self.returnCode != 0:
            print ("[PIPE:] STOP! Job '%s' caused non-zero exit code (%s)" % (self.cmd, self.returnCode))
        elif cleanup: ## only upon successful completion
            self.cleanup()
        if 'process' in vars(self):
            del self.process
        return self.returnCode==0
    def mark4cleanup(self, tmpfile, condition=None):
        self.marked4deletion.append((tmpfile, condition))
    def cleanup(self):
        for tmpfile, condition in self.marked4deletion: #usually just one file/job, typically the job input
            if condition is None or (os.path.exists(condition) and os.stat(condition).st_size > 100000):
                subprocess.Popen('rm %s'%tmpfile, shell=True) # can handle patterns!
                #os.remove(tmpfile)
def pptime(s):
    hours, remainder = divmod(s, 3600)
    minutes, seconds = divmod(remainder, 60)
    return '%02d:%02d:%02d' % (hours, minutes, seconds)

if __name__ == "__main__":
    restartFromCheckpoint = False
    ## set to False, if you want to start from the beginning

    if len(sys.argv[-1]) < 3: ## interprete small numbers as job ids
        sampleNr = int(sys.argv[-1])
        pipe = Pipeline(sampleNr=sampleNr, refShort='hg19', singleLane=False) 
    else: ## interprete 5 or 6digit numbers as sample ids
        sample = sys.argv[-1]
        pipe = Pipeline(sample=sample, refShort='hg19', singleLane=False)
        
    lastSuccJobIdx = pipelineOrder.index('ApplyBQSR')
    if True:
        lastSuccJobIdx = pipe.getLastCheckpoint() if restartFromCheckpoint else -1

        if pipelineOrder.index('trimming') > lastSuccJobIdx:
            pipe.trim()  #pattern="*L002*.fastq.gz") ## example for just trimming lane2 files
            pipe.report()
            pipe.fastQC(pipe.dirs['trimOutdir'], pipe.dirs['fastqOutdir'])

        ## BWA
        if pipelineOrder.index('bwa') > lastSuccJobIdx:
            pipe.align()
            pipe.report()

        ## PICARD
        if pipelineOrder.index('SortSam') > lastSuccJobIdx:
            pipe.sortSam(jvmSpaceDivide=False) 
            pipe.report()

        if pipelineOrder.index('MergeBam') > lastSuccJobIdx:
            pipe.mergeBam()
            pipe.report()

        if pipelineOrder.index('MarkDuplicates') > lastSuccJobIdx:
            pipe.markDuplicates() #126G
            pipe.report()
            pipe.qualimap()

        if pipelineOrder.index('BuildBamIndex') > lastSuccJobIdx:
            pipe.buildBamIndex()
            pipe.report()
    if True:
        ## GATK
        if pipelineOrder.index('BaseRecalibrator') > lastSuccJobIdx:
            pipe.baseRecalibrator()
            pipe.report()

        if pipelineOrder.index('ApplyBQSR') > lastSuccJobIdx:
            pipe.applyBQSR()
            pipe.report()

        ##
        if pipelineOrder.index('HaplotypeCaller') > lastSuccJobIdx:
            pipe.haplotypeCaller()
            pipe.report()
            
    if False:
        ## Before variant recal: need to do combineGVCFs, GenotypeGVCFs -> VCF
        pipe.filetype='gvcf'
        pipe.zipped = '.gz'
        pipe.variantRecalibrator() ## only works on proper VCF!!!
        pipe.report()
    
    #pipe.applyVQSR()
    #pipe.variantRecalibrator(pattern='*_recal_SNP', mode='INDEL')
    #pipe.applyVQSR(pattern='*_recal_SNP', mode='INDEL')


