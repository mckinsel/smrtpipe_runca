"""


        a congeries of iridescent globes



"""

import logging
import os

from SMRTpipe.modules.P_Module import P_Module
from SMRTpipe.engine.DistributableTasks import DistributableTask
from SMRTpipe.engine.SmrtPipeWorkflow import SmrtWorkflowError
from SMRTpipe.engine.SmrtPipeTasks import task
from SMRTpipe.cluster.Scatter import ScatterFofn
from SMRTpipe.cluster.Gather import SentinelGather, GatherWithCat
from SMRTpipe.engine.common import USE_GLOBAL_NPROC
from SMRTpipe.engine.SmrtPipeFiles import SMRTDataFile, SMRTFile
from SMRTpipe.modules.P_CeleraAssemblerBase import P_CeleraAssemblerBase
from P_PreAssemblerDagcon import dconCorrectedFastq, dconCorrectedFasta

log = logging.getLogger(__name__)

################################################################################
# The location of the Celera Assembler binaries needs to be set as
# an environment variable for now, since they're not kept in the
# usual place in the build.
try:
    CA_BIN_DIR = os.environ['CELERA_ASSEMBLER_DIR']
except KeyError as e:
    print 'Must define the env variable CELERA_ASSEMBLER_DIR'
    raise(e)
################################################################################

################################################################################
################################################################################
# Files produced by Celera Assembler
################################################################################
################################################################################

# The frg file that's used at the input to Celera Assembler. It mostly just
# points to the fastq file produced by pbdagcon
corrReadsFrg = SMRTFile("data/corrected.frg")

# Celera's 'gatekeeper store'. This maintains the reads themselves as well as
# the substring of the reads that have passed quality filtering. It's edited
# by many Celera processes, so you can't really rerun CA without creating the
# gkpStore anew.
gkpStore  = SMRTFile("data/ca_gkp_store")

# A summary of the initial, base quality trimming. Nobody cares about this
# file.
initialTrimSummary = SMRTFile("data/ca_initial_trim.summary")
# The log file for initial trimming. I don't think it's important either, but
# here we are.
initialTrimLog = SMRTFile("data/ca_initial_trim.log")

# Meryl outputs two files that presumably contain kmer count information about
# the input fastq file. Meryl assigns the extentions, it just wants the
# basenames.
merylDatOutput = SMRTFile("data/ca_meryl_output.mcdat")
merylIdxOutput = SMRTFile("data/ca_meryl_output.mcidx")

# Meryl will get a threshold for overrepresented (?) kmers. This gets written
# to this file, which just contains an integer and that's it.
merylThresholdFile = SMRTFile("data/ca_meryl_mer_thresh.out")
# Finally, Meryl will write the kmers themselves in this fasta file.
# Subsequent steps in the CA pipeline use this file.
merylNmersFasta = SMRTFile("data/ca_meryl_nmers.fasta")

# CA has a binary called overlapInCore which will align the reads against
# themselves, recording the overlaps. This task is split by aligning subsets of
# the reads against other subsets. The ovl_partitions file describes how they
# are to be divided
ovlPartitionsFile = SMRTFile("data/ca_overlap_partitions.txt")

obtOvlDir = SMRTFile("data/ca_obt_overlaps")

# We need a file to track when all the overlap jobs are done. This is it! 
obtOvlSentinel = SMRTFile("data/.ca_obt_ovl_sentinel")

# The list of all the overlap results files
obtOvlList = SMRTFile("data/ca_obt_ovl.fofn")

# The Celera overlapStore
obtOvlStore = SMRTFile("data/ca_obtOvlStore")

# And so on and so forth...
finalTrimSummary = SMRTFile("data/ca_final_trim.summary")
chimeraSummary = SMRTFile("data/ca_chimera.summary")

ovlOvlDir = SMRTFile("data/ca_ovl_overlaps")
ovlOvlSentinel = SMRTFile("data/.ca_ovl_ovl_sentinel")
ovlOvlList = SMRTFile("data/ca_ovl_ovl.fofn")
ovlOvlStore = SMRTFile("data/ca_ovlOvlStore")

fragCorr = SMRTFile("data/ca_frag_corr")
fragCorrList = SMRTFile("data/ca_frg_corr.fofn")
catFragCorr = SMRTFile("data/ca_cat_frag_corr")
erateFile = SMRTFile("data/ca_erate")
erateList = SMRTFile("data/ca_erate.fofn")
catErate = SMRTFile("data/ca_cat_erate")
ovlUpdateSentinel = SMRTFile("data/.ca_ovl_update_sentinel")
tigStore = SMRTFile("data/ca_tig_store")
unitigList = SMRTFile("data/ca_utg.list")
utgConsensus = SMRTFile("data/draft_assembly.fasta")

################################################################################
################################################################################
# Celera Assembler commands
################################################################################
################################################################################

# Create the gatekeeper store
cmd_create_gkp = "{ca_bin}/gatekeeper -o {gkp}  -T -F {frg}"

# Perform the initial quality-based trimming
cmd_initial_trim = "{ca_bin}/initialTrim -log {log} -frg {gkp} > {summary}"

# Several commands to get the meryl nmers
cmd_meryl_first = ("{ca_bin}/meryl -B -C -v -m {mer_size} -memory {memory_size} -threads {threads} "
                   "-c {mer_comp} -L 2 -s {gkp}:chain -o {mer_out}")
cmd_mer_threshold= "{ca_bin}/estimate-mer-threshold -m {mer_out} > {mer_thresh_out}"
cmd_meryl_second = "{ca_bin}/meryl -Dt -n $(cat {mer_thresh}) -s {mer_out} > {nmers_fasta}"

# Get the partitioning for the overlap jobs
cmd_ovl_partition = 'calculate_ovl_ranges.py $(grep -c "^>" {corrected_fasta}) {max_chunks} > {ovl_partitions}'

# Do the actual overlapping
cmd_overlap = ('run_overlap_in_core.py --ovl_partitions {ovl_partitions} '
                   '--hashbits {hash_bits} --gkp_store {gkp} '
                   '--nmers_fasta {nmers} --ca_bin_dir {ca_bin_dir} '
                   '--output_dir {out_dir}')

# Create overlapStore, which contains the results of the overlapping
cmd_build_overlap_store = "{ca_bin}/overlapStoreBuild -o {ovl_store} -g {gkp} -M {max_mem} -L {ovl_list}" 

# Perform final trimming, this looks at overlap patterns
cmd_final_trim = "{ca_bin}/finalTrim -G {gkp} -O {ovl_store} -e 0.03 -E 3.25 -o {final_trim_summary}"

# Run the chimera and spur detection steps
cmd_chimera = "{ca_bin}/chimera -G {gkp} -O {ovl_store} -e 0.06 -E 1.1 -o {chimera_summary}"

# Some mysterious correction steps
cmd_frg_correct = '{ca_bin}/correct-frags -t 2 -S {ovl_store} -o {frg_corr_out} {gkp} 1 $(grep -c "^>" {corrected_fasta})'
cmd_cat_corrects = "{ca_bin}/cat-corrects -L {frg_corr_list} -o {cat_frg_corr}"
cmd_correct_olaps = '{ca_bin}/correct-olaps -S {ovl_store} -e {erate} {gkp} {frg_corr} 1 $(grep -c "^>" {corrected_fasta})'
cmd_cat_erate = '{ca_bin}/cat-erates -L {erate_list} -o {cat_erates}'
cmd_update_ovl_store = '{ca_bin}/overlapStore -u {ovl_store} {cat_erates}'

# Run Celera's unitigger
cmd_bogart = ('{ca_bin}/bogart -O {ovl_store} -G {gkp} -T {tig_store} -B 100000 '
              '-eg 0.03  -Eg 3.25  -em 0.045  -Em 5.25 -o {bogart_root}')

# Create the tigStore from which we'll get the draft assembly
cmd_get_unitigs = ('{ca_bin}/tigStore -g {gkp} -t {tig_store} 1 -d properties -U | '
                   'awk \'BEGIN{{t=0}}$1=="numFrags"{{if($2>1){{print t, $2}}t++}}\' | '
                   'sort -nrk2,2 > {unitig_list}')

class P_RunCA(P_Module):
    
    def validateSettings( self ):
        errors = P_Module.validateSettings(self)

        self.merylMerSize = self.setting('meryl_mer_size', '14')
        self.merylMemorySize = self.setting('meryl_memory_size', '1374') 
        self.merylCompression = self.setting('meryl_compression', '0')
        
        self.ovlHashBits = self.setting("ovl_hash_bits", '22')

        return errors
    
    @task(inputs= {'corrReadsFastq': dconCorrectedFastq},
          outputs={'corrReadsFrg': corrReadsFrg})
    def genFrgFile(self, files):
        ifile = files.corrReadsFastq.path
        ofile = files.corrReadsFrg.path
        cmd = 'fastqToCA -technology sanger -type sanger '
        cmd += '-reads %s -libraryname %s > %s' \
            % (ifile, "pacbio", ofile)
        return cmd


    @task(inputs= {'corrReadsFrg': corrReadsFrg},
          outputs={'gkpStore': gkpStore},
          nproc=1)
    def createGkpStore(self, files):
        cmd = cmd_create_gkp.format(ca_bin=CA_BIN_DIR, gkp=files.gkpStore.path,
                                    frg=files.corrReadsFrg.path)
        return cmd 

    @task(inputs= {'gkpStore': gkpStore},
          outputs={'initialTrimSummary': initialTrimSummary,
                   'initialTrimLog': initialTrimLog},
          nproc=1)
    def initialTrim(self, files):
        cmd = cmd_initial_trim.format(ca_bin=CA_BIN_DIR, gkp=files.gkpStore.path,
                                      log=files.initialTrimLog.path,
                                      summary=files.initialTrimSummary.path)
        return cmd

    @task(inputs= {'gkpStore': gkpStore},
          outputs={'merylDatOutput': merylDatOutput,
                   'merylIdxOutput': merylIdxOutput,
                   'merylThresholdFile': merylThresholdFile,
                   'merylNmersFasta': merylNmersFasta},
          nproc=USE_GLOBAL_NPROC)
    def meryl(self, files):

        meryl_output_root, ext = os.path.splitext(files.merylDatOutput.path)
        cmds = []

        meryl_cmd = cmd_meryl_first.format(ca_bin=CA_BIN_DIR, mer_size=self.merylMerSize,
                                           memory_size=self.merylMemorySize, threads=files.nproc,
                                           mer_comp=self.merylCompression, gkp=files.gkpStore.path,
                                           mer_out=meryl_output_root)
        cmds.append(meryl_cmd)

        thresh_cmd = cmd_mer_threshold.format(ca_bin=CA_BIN_DIR, mer_out=meryl_output_root,
                                              mer_thresh_out=files.merylThresholdFile.path)
        cmds.append(thresh_cmd)

        meryl_nmers_cmd = cmd_meryl_second.format(
            ca_bin=CA_BIN_DIR, mer_thresh=files.merylThresholdFile.path,
            mer_out=meryl_output_root, nmers_fasta=files.merylNmersFasta.path)
        cmds.append(meryl_nmers_cmd)

        return cmds
    
    @task(inputs= {'dconCorrectedFasta': dconCorrectedFasta},
          outputs={'ovlPartitionsFile': ovlPartitionsFile},
          nproc=1,
          localOnly=True)
    def overlap_partition(self, files):
        out_dir = os.path.dirname(files.ovlPartitionsFile.path)

        cmd = cmd_ovl_partition.format(corrected_fasta=files.dconCorrectedFasta.path,
                                       max_chunks=15,
                                       ovl_partitions=files.ovlPartitionsFile.path)
        return cmd

    @task(inputs= {'gkpStore': gkpStore,
                   'ovlPartitionsFile': ovlPartitionsFile,
                   'merylNmersFasta': merylNmersFasta,
                   'obtOvlDir': obtOvlDir},
          outputs={'obtOvlSentinel': obtOvlSentinel},
          taskType=DistributableTask,
          chunkFunc=lambda x: 15,
          scatters=[ScatterFofn(ovlPartitionsFile)],
          gathers=[SentinelGather(obtOvlSentinel)])
    def obtOverlap(self, files):
        if not os.path.exists(files.obtOvlDir.path):
            os.mkdir(files.obtOvlDir.path)
        cmds = [] 
        cmd = cmd_overlap.format(ovl_partitions=files.ovlPartitionsFile.path, 
                                 gkp=files.gkpStore.path, nmers=files.merylNmersFasta.path,
                                 ca_bin_dir=CA_BIN_DIR, hash_bits=self.ovlHashBits,
                                 out_dir=files.obtOvlDir.path)
        cmd += ' --do_partial'
        cmds.append(cmd)
        cmds.append("touch {o}".format(o=files.obtOvlSentinel.path))
        return cmds

    @task(inputs= {'obtOvlSentinel': obtOvlSentinel,
                   'obtOvlDir': obtOvlDir},
          outputs={'obtOvlList': obtOvlList},
          nproc=1,
          localOnly=True)
    def gatherObtOvls(self, files):
        data_dir = self._context.getDataService().getDataPath()
        cmd = "for i in $(ls {d}/*.ovb.gz); do readlink -f $i >> {c}; done".format(
            d=files.obtOvlDir.path, c=files.obtOvlList.path)
        
        return cmd

    @task(inputs= {'gkpStore': gkpStore,
                   'obtOvlList': obtOvlList},
          outputs={'obtOvlStore': obtOvlStore},
          nproc=1)
    def buildObtOverlapStore(self, files):
        cmd = cmd_build_overlap_store.format(ca_bin=CA_BIN_DIR, gkp=files.gkpStore.path,
                                             max_mem=self.merylMemorySize, ovl_list=files.obtOvlList.path,
                                             ovl_store=files.obtOvlStore.path) 
        cmd += ' -obt'
        return cmd

    @task(inputs= {'gkpStore': gkpStore,
                   'obtOvlStore': obtOvlStore},
          outputs={'finalTrimSummary': finalTrimSummary})
    def finalTrim(self, files):
        final_trim_root, ext = os.path.splitext(files.finalTrimSummary.path)
        cmd = cmd_final_trim.format(ca_bin=CA_BIN_DIR, gkp=files.gkpStore.path,
                                    ovl_store=files.obtOvlStore.path,
                                    final_trim_summary=final_trim_root)
        return cmd


    @task(inputs= {'gkpStore': gkpStore,
                   'obtOvlStore': obtOvlStore,
                   'finalTrimSummary': finalTrimSummary},
          outputs={'chimeraSummary': chimeraSummary})
    def chimera(self, files):
        chimera_root, ext = os.path.splitext(files.chimeraSummary.path)
        cmd = cmd_chimera.format(ca_bin=CA_BIN_DIR, gkp=files.gkpStore.path,
                                 ovl_store=files.obtOvlStore.path,
                                 chimera_summary=chimera_root)
        return cmd

    @task(inputs= {'gkpStore': gkpStore,
                   'ovlPartitionsFile': ovlPartitionsFile,
                   'merylNmersFasta': merylNmersFasta,
                   'ovlOvlDir': ovlOvlDir,
                   'chimeraSummary': chimeraSummary},
          outputs={'ovlOvlSentinel': ovlOvlSentinel},
          taskType=DistributableTask,
          chunkFunc=lambda x: 15,
          scatters=[ScatterFofn(ovlPartitionsFile)],
          gathers=[SentinelGather(ovlOvlSentinel)])
    def ovlOverlap(self, files):
        if not os.path.exists(files.ovlOvlDir.path):
            os.mkdir(files.ovlOvlDir.path)
        cmds = [] 
        cmd = cmd_overlap.format(ovl_partitions=files.ovlPartitionsFile.path, 
                                 gkp=files.gkpStore.path, nmers=files.merylNmersFasta.path,
                                 ca_bin_dir=CA_BIN_DIR, hash_bits=self.ovlHashBits,
                                 out_dir=files.ovlOvlDir.path)
        cmds.append(cmd)
        cmds.append("touch {o}".format(o=files.ovlOvlSentinel.path))
        return cmds

    @task(inputs= {'ovlOvlSentinel': ovlOvlSentinel,
                   'ovlOvlDir': ovlOvlDir},
          outputs={'ovlOvlList': ovlOvlList},
          nproc=1,
          localOnly=True)
    def gatherOvlOvls(self, files):
        data_dir = self._context.getDataService().getDataPath()
        cmd = "for i in $(ls {d}/*.ovb.gz); do readlink -f $i >> {c}; done".format(
            d=files.ovlOvlDir.path, c=files.ovlOvlList.path)
        
        return cmd

    @task(inputs= {'gkpStore': gkpStore,
                   'ovlOvlList': ovlOvlList},
          outputs={'ovlOvlStore': ovlOvlStore},
          nproc=1)
    def buildOvlOverlapStore(self, files):
        cmd = cmd_build_overlap_store.format(ca_bin=CA_BIN_DIR, gkp=files.gkpStore.path,
                                             max_mem=self.merylMemorySize, ovl_list=files.ovlOvlList.path,
                                             ovl_store=files.ovlOvlStore.path) 
        return cmd

    @task(inputs= {'ovlOvlStore': ovlOvlStore,
                   'gkpStore': gkpStore,
                   'dconCorrectedFasta': dconCorrectedFasta},
          outputs={'fragCorr': fragCorr,
                   'fragCorrList': fragCorrList,
                   'catFragCorr': catFragCorr})
    def fragCorrect(self, files):
        cmds = []
        cmd = cmd_frg_correct.format(ca_bin=CA_BIN_DIR, gkp=files.gkpStore.path,
                                     ovl_store=files.ovlOvlStore.path,
                                     corrected_fasta=files.dconCorrectedFasta.path,
                                     frg_corr_out=files.fragCorr.path)
        cmds.append(cmd)

        list_cmd = "echo {f} > {l}".format(f=files.fragCorr.path, l=files.fragCorrList.path)
        cmds.append(list_cmd)
        cat_cmd = cmd_cat_corrects.format(ca_bin=CA_BIN_DIR, frg_corr_list=files.fragCorrList.path,
                                          cat_frg_corr=files.catFragCorr.path)
        cmds.append(cat_cmd)
        return cmds


    @task(inputs= {'ovlOvlStore': ovlOvlStore,
                   'gkpStore': gkpStore,
                   'catFragCorr': catFragCorr,
                   'dconCorrectedFasta': dconCorrectedFasta},
          outputs={'erateFile': erateFile,
                   'erateList': erateList,
                   'catErate': catErate})
    def olapCorrect(self, files):
        cmds = []
        cmd = cmd_correct_olaps.format(ca_bin=CA_BIN_DIR, ovl_store=files.ovlOvlStore.path,
                                       gkp=files.gkpStore.path, frg_corr=files.catFragCorr.path,
                                       corrected_fasta=files.dconCorrectedFasta.path,
                                       erate=files.erateFile.path)
        cmds.append(cmd)
        list_cmd = "echo {f} > {l}".format(f=files.erateFile.path, l=files.erateList.path)
        cmds.append(list_cmd)
        
        cat_cmd = cmd_cat_erate.format(ca_bin=CA_BIN_DIR, erate_list=files.erateList.path,
                                       cat_erates=files.catErate.path)
        cmds.append(cat_cmd)
        return cmds

    @task(inputs= {'ovlOvlStore': ovlOvlStore,
                   'catErate': catErate},
          outputs={'ovlUpdateSentinel': ovlUpdateSentinel})
    def updateOvlStore(self, files):
        cmds = []
        cmd = cmd_update_ovl_store.format(ca_bin=CA_BIN_DIR, ovl_store=files.ovlOvlStore.path,
                                          cat_erates=files.catErate.path)
        cmds.append(cmd)
        cmds.append("touch {o}".format(o=files.ovlUpdateSentinel.path))
        return cmds

    @task(inputs= {'ovlUpdateSentinel': ovlUpdateSentinel,
                   'gkpStore': gkpStore,
                   'ovlOvlStore': ovlOvlStore},
          outputs={'tigStore': tigStore})
    def bogart(self, files):

        data_dir = self._context.getDataService().getDataPath()
        bogart_root = os.path.join(data_dir, "ca_bogart")

        cmd = cmd_bogart.format(ca_bin=CA_BIN_DIR, gkp=files.gkpStore.path,
                                ovl_store=files.ovlOvlStore.path,
                                tig_store=files.tigStore.path,
                                bogart_root=bogart_root)
        return cmd

    @task(inputs= {'tigStore': tigStore,
                   'gkpStore': gkpStore},
          outputs={'unitigList': unitigList})
    def getUnitigs(self, files):
        
        cmd = cmd_get_unitigs.format(ca_bin=CA_BIN_DIR, gkp=files.gkpStore.path,
                                     tig_store=files.tigStore.path,
                                     unitig_list=files.unitigList.path)
        return cmd
    
    @task(inputs={'unitigList': unitigList,
                  'tigStore': tigStore,
                  'gkpStore': gkpStore},
          outputs={'utgConsensus': utgConsensus},
          nproc=USE_GLOBAL_NPROC,
          taskType=DistributableTask,
          scatters=[ScatterFofn(unitigList)],
          chunkFunc=lambda x: 5,
          gathers=[GatherWithCat(utgConsensus)],
          localOnly=False)
    def unitigConsensus(self, files):
        data_dir = self._context.getDataService().getDataPath()
        ca_prefix = os.path.join(data_dir, 'ca_utgcon')
        cmds = []
        cmds.append("tmp=$(TMPDIR=%s mktemp -d)" % self._context.config['TMP'])
        cmds.append("tmp=$tmp tig=%s gkp=%s utg=%s nproc=%s cns=%s pbutgcns_wf.sh"
            % (files.tigStore.path, files.gkpStore.path, files.unitigList.path,
               files.nproc, files.utgConsensus.path))
        return cmds
