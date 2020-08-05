
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
    outdir = "/research/btc_bioinformatic/operations/scratch/"
    tmpdir = "/research/btc_bioinformatic/operations/scratch/tmp/"
    
    applicationDir = '/research/btc_bioinformatic/operations/'
    rawData = '/research/btc_bioinformatic/results/'
    resourceDir = '/research/btc_bioinformatic/operations/data_masdar/'
    Mills = "%s/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf" % resourceDir
    G1000_indels = "%s/1000G_phase1.indels.hg19.sites.vcf"% resourceDir
    G1000_snps = "%s/1000G_phase1.snps.high_confidence.hg19.sites.vcf"% resourceDir
    Hapmap = "%s/hapmap_3.3.hg19.sites.vcf"% resourceDir
    Omni = "%s/1000G_omni2.5.hg19.sites.vcf"% resourceDir
    dbsnp = "%s/dbsnp_138.hg19.vcf"% resourceDir  # could change to a most recent dbsnp
    TR = "%s/truseq-exome-targeted-regions-manifest-v1-2a.bed"% resourceDir  # Target region for sequencing WES
    jobOutputDir = '%s/Pipeline/JobOutputs' % applicationDir


    rawDataDir = '%s/WES_NF' % rawData

    refDirLookup = {'hg19':'%s/data_masdar/ucsc.hg19.fasta' % applicationDir, 
                    'hg38':'%s/AuxData/hg38/Homo_sapiens_assembly38.fasta' % applicationDir}


class Pipeline2:
    def __init__(self, refShort='hg19', vcfdir='', base='/research/btc_bioinformatic/results/gathered_vcfs'):
        self.base = base
        self.gatheredvcf = base + '.vcf.gz' 
        self.outvcf = {'SNP': base + '_recal_snp.vcf.gz',
                       'INDEL': base + '_recal_snp_indel.vcf.gz'}
        self.refShort = refShort
        self.vcfdir = vcfdir
        self.ref=refDirLookup[refShort]
        self.completedJobs = []
        self.tranches = '_recal_%s.tranches'
        self.recal = '_recal_%s.recal'
        
    def variantRecalibrator(self, jvmSpace=32, mode='SNP'):

        recal_output = self.base + self.recal % (mode)
        tranches_output = self.base + self.tranches % mode
        plots_output = self.base + '_recal_%s.rplots.R' % mode
        G1000 = {'SNP': G1000_snps, 'INDEL': G1000_indels}[mode]
        if mode=='SNP':
            vcfFile = self.gatheredvcf
            modeSpecificResources = ['--resource hapmap,known=false,training=true,truth=true,prior=15.0:%s' % Hapmap,
                                     '--resource omni,known=false,training=true,truth=true,prior=12.0:%s' % Omni]
        elif mode=='INDEL':
            vcfFile = self.outvcf['SNP']
            modeSpecificResources = ['--resource mills,known=false,training=true,truth=true,prior=12.0:%s' % Mills]

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

    def applyVQSR(self, jvmSpace=32, filterLevel='99.0', mode='SNP'):
        if mode=='SNP':
            vcfFile = self.gatheredvcf
        elif mode=='INDEL':
            vcfFile = self.outvcf['SNP']
        output = self.outvcf[mode]
        recal_input = self.base + self.recal % mode
        tranches_input = self.base + self.tranches % mode 
        recCmd = ['gatk', '--java-options', '"-Xmx%sG"' % jvmSpace, 'ApplyVQSR',
                  '-R', self.ref,
                  '-V', vcfFile,
                  '-O', output,
                  '-ts-filter-level', filterLevel, # --ts_filter_level
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
    import sys
    regionID = int(sys.argv[-1])-1 if len(sys.argv)>1 else 0
    ## last argument is taken as regionID, should be a number between 1 and 322 (or however many vcfs there are in the vcfdir)
    pipe = Pipeline2(refShort='hg19', vcfdir='/research/btc_bioinformatic/results/GenotypeGVCFs')

    pipe.variantRecalibrator(mode='SNP') ## only works on proper VCF!!!
    #pipe.report()
    pipe.applyVQSR(mode='SNP')
    pipe.variantRecalibrator(mode='INDEL') ## only works on proper VCF!!!
    #pipe.report()
    pipe.applyVQSR(mode='INDEL')    
    

