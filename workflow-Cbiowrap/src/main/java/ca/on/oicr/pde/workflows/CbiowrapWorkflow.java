package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import org.apache.commons.io.FilenameUtils;
import org.mortbay.log.Log;

/**
 * <p>
 * For more information on developing workflows, see the documentation at
 * <a href="http://seqware.github.io/docs/6-pipeline/java-workflows/">SeqWare
 * Java Workflows</a>.</p>
 *
 * Quick reference for the order of methods called: 1. setupDirectory 2.
 * setupFiles 3. setupWorkflow 4. setupEnvironment 5. buildWorkflow
 *
 * See the SeqWare API for
 * <a href="http://seqware.github.io/javadoc/stable/apidocs/net/sourceforge/seqware/pipeline/workflowV2/AbstractWorkflowDataModel.html#setupDirectory%28%29">AbstractWorkflowDataModel</a>
 * for more information.
 */
public class CbiowrapWorkflow extends OicrWorkflow {

    //dir
    private String dataDir, tmpDir;
    private String outDir;
    private Set<String> keySet;
    // Input Data
    private String inputCommaSeparatedMafs;
    private ArrayList<String> inputMafs = new ArrayList<String> ();
    private int nMafs;
    private String inputCommaSeparatedCopyNumberSegs;
    private ArrayList<String> inputSegs = new ArrayList<String> ();
    private int nSegs;
    private String inputCommaSeparatedRSEMCounts;
    private ArrayList<String> inputRcounts = new ArrayList<String> ();
    private int nRCs;
    private String inputCommaSeparatedSTARrtabs;
    private ArrayList<String> inputRtabs = new ArrayList<String> ();

    
    
    // cbiwrap Rscript
    private String cbioWrapper;
    private String rPath;
    private String cbioWrapPath;
    
    // reference files
    private String hotspotGenesFile = "/.mounts/labs/TGL/gsi/databases/20161121_Chang_hotspot_list.tsv";
    private String oncoKBFile = "/.mounts/labs/TGL/gsi/databases/20170412_oncoKB.tsv";
    private String ensembleConversonFile = "/.mounts/labs/TGL/gsi/databases/ensemble_conversion.txt";
    private String blackList ="/.mounts/labs/TGL/gsi/databases/blacklist.txt";
    
    // additional references (for future use)
    private String additionalReferences = null;
    
    // optional CSV file for human id mapping
    private String mastersheetCSV = null;
    
    
    private boolean manualOutput;
    private String queue;
    private int cbiowrapMem;
    
    //

    
    // metatypes
    private final String TXT_GZ_METATYPE="application/text-gz";
    private final String FOLDER_GZ_METATYPE="application/tar-gz";
    private final String TEXT_METATYPE="text/plain";
    
    

    private void init() {
        try {
            //dir
            dataDir = "data/";
            tmpDir = getProperty("tmp_dir");
            
            // read scripts
            rPath = getProperty("rpath");
            cbioWrapPath = getProperty("cbiowrap_basedir");
            cbioWrapper = cbioWrapPath + "/R/wrapper.r";

            /// read all params
            inputCommaSeparatedMafs = getProperty("input_mafs");
            inputCommaSeparatedCopyNumberSegs = getProperty("input_segs");
            inputCommaSeparatedRSEMCounts = getOptionalProperty("input_rsem_counts", "blank");
            inputCommaSeparatedSTARrtabs = getOptionalProperty("input_star_rtabs", "blank");

            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            cbiowrapMem = Integer.parseInt(getProperty("cbiowrap_mem"));

        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void setupDirectory() {
        init();
        this.addDirectory(dataDir);
        this.addDirectory(tmpDir);
        if (!dataDir.endsWith("/")) {
            dataDir += "/";
        }
        if (!tmpDir.endsWith("/")) {
            tmpDir += "/";
        }
    }

    @Override
    public Map<String, SqwFile> setupFiles() {
        //provision in MAFs
        String[] allMAfs = inputCommaSeparatedMafs.split(",");
        nMafs = allMAfs.length;
        for (int imx = 0; imx < nMafs; imx++){
            String iMaf = allMAfs[imx];
            SqwFile file0 = this.createFile("inputMAFs" + Integer.toString(imx));
            file0.setSourcePath(iMaf);
            file0.setType(TXT_GZ_METATYPE);
            file0.setIsInput(true);  
        }
        //provision in SEGs
        String[] allSegs = inputCommaSeparatedCopyNumberSegs.split(",");
        nSegs = allSegs.length;
        for (int isx = 0; isx < nSegs; isx++){
            String iSeg = allSegs[isx];
            SqwFile file1 = this.createFile("inputSEGs" + Integer.toString(isx));
            file1.setSourcePath(iSeg);
            file1.setType(TEXT_METATYPE);
            file1.setIsInput(true);  
        }
         //
        if (!this.inputCommaSeparatedRSEMCounts.equals("blank") || !this.inputCommaSeparatedSTARrtabs.equals("blank")) {
            Map<String, List<String>> inputFileMap = this.getRsemStarMap(this.inputCommaSeparatedRSEMCounts, this.inputCommaSeparatedSTARrtabs);
            this.keySet = inputFileMap.keySet();
            for (String key : keySet) {
                SqwFile file0 = this.createFile("RSEM_" + key);
                SqwFile file1 = this.createFile("STAR_" + key);
                String rsemFile = inputFileMap.get(key).get(0);
                String starFile = inputFileMap.get(key).get(1);
                file0.setSourcePath(rsemFile);
                file0.setType(TEXT_METATYPE);
                file0.setIsInput(true);
                file1.setSourcePath(starFile);
                file1.setType(TEXT_METATYPE);
                file1.setIsInput(true);
            }
        }
        
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {
        Job parentJob = null;
        
        Date date = Calendar.getInstance().getTime();
        DateFormat dateFormat = new SimpleDateFormat("yyyymmdd");
        String outputFileNamePrefix = dateFormat.format(date);
        // read MAFs
        for (int imx = 0; imx <= nMafs; imx ++){
            String inMAF = getFiles().get("inputMAFs" + Integer.toString(imx)).getProvisionedPath();
            inputMafs.add(inMAF);
        }
        // post process MAFs (copy to tmp directory)
        String mafDir = this.dataDir + "MAF";
        String inMAFs = inputMafs.toString().replace(", ", ",").replaceAll("[\\[.\\]]", "");
        String combinedMaf = this.tmpDir + outputFileNamePrefix + "_combined_raw.maf.txt";
        Job createMafDir = linkMAFs(mafDir, inMAFs, combinedMaf);
        parentJob = createMafDir;
        
        // read Segs
        for (int isx = 0; isx <= nSegs; isx++){
            String inSEG = getFiles().get("inputSEGs" + Integer.toString(isx)).getProvisionedPath();
            inputSegs.add(inSEG);
        }
        // post process Segs into a single file
        String combinedSegFile = this.tmpDir + outputFileNamePrefix  + "_copynumber.seg";
        String inSegFiles = inputSegs.toString().replace(", ", ",").replaceAll("[\\[.\\]]", "");
        Job generateCombinedCopySeg = copySegProcess(inSegFiles, combinedSegFile);
        // provision out combined seg file
        SqwFile copySegSqw = createOutputFile(combinedSegFile, TEXT_METATYPE, this.manualOutput);
        copySegSqw.getAnnotations().put("SEG", "cbiowrap");
        generateCombinedCopySeg.addFile(copySegSqw);
        
        generateCombinedCopySeg.addParent(parentJob);
        parentJob = generateCombinedCopySeg;
        
        String postprocRC = null;
        if (!this.inputCommaSeparatedRSEMCounts.equals("blank") || !this.inputCommaSeparatedSTARrtabs.equals("blank")) {
            // read RTABs
            for (int itx = 0; itx <= nRCs; itx++) {
                String inRT = getFiles().get("inputRTABs" + Integer.toString(itx)).getProvisionedPath();
                inputRtabs.add(inRT);
            }
            String csInRTs = inputRtabs.toString().replace(", ", ",").replaceAll("[\\[.\\]]", "");
            
            // read Rcounts
            for (int irx = 0; irx <= nRCs; irx++) {
                String inRC = getFiles().get("inputRCOUNTs" + Integer.toString(irx)).getProvisionedPath();
                inputRcounts.add(inRC);
            }
            // post process Rcounts 
            String csInRCs = inputRcounts.toString().replace(", ", ",").replaceAll("[\\[.\\]]", "");
            postprocRC = this.dataDir + outputFileNamePrefix + "_RCOUNT.txt";
            Job postProcRSEM = postProcessRSEM(csInRCs, csInRTs, postprocRC);
            // provision out combined Rcounts file
            SqwFile rcountsSqw = createOutputFile(postprocRC, TEXT_METATYPE, this.manualOutput);
            rcountsSqw.getAnnotations().put("RCOUNTS", "cbiowrap");
            postProcRSEM.addFile(rcountsSqw);
            // provision out combine tpm file
            String tpm = postprocRC.replace("_RCOUNT.txt", "_GENES_TPM.txt");
            SqwFile tpmSqw = createOutputFile(tpm, TEXT_METATYPE, this.manualOutput);
            tpmSqw.getAnnotations().put("TPM", "cbiowrap");
            postProcRSEM.addFile(tpmSqw);
            // provision out combine fpkm file
            String fpkm = postprocRC.replace("_RCOUNT.txt", "_GENES_FPKM.txt");
            SqwFile fpkmSqw = createOutputFile(fpkm, TEXT_METATYPE, this.manualOutput);
            fpkmSqw.getAnnotations().put("FPKM", "cbiowrap");
            postProcRSEM.addFile(fpkmSqw);
            // provision out combine count file
            String count = postprocRC.replace("_RCOUNT.txt", "_GENES_COUNT.txt");
            SqwFile countSqw = createOutputFile(count, TEXT_METATYPE, this.manualOutput);
            countSqw.getAnnotations().put("COUNT", "cbiowrap");
            postProcRSEM.addFile(countSqw);
            
            postProcRSEM.addParent(parentJob);
            parentJob = postProcRSEM;
        }
        
        // run cbiowrap
        String cbioWrapDir = this.dataDir + "wrapper_" + outputFileNamePrefix;
        Job launchCBioWrap = runCbioWrap(combinedMaf, combinedSegFile, postprocRC, cbioWrapDir);
        launchCBioWrap.addParent(parentJob);
        parentJob = launchCBioWrap;
        // provison out cbiowrap dir
        String cbioWrapDirGz = cbioWrapDir + ".tar.gz";
        SqwFile cbioWrapDirGzSqw = createOutputFile(cbioWrapDirGz, FOLDER_GZ_METATYPE, this.manualOutput);
        cbioWrapDirGzSqw.getAnnotations().put("cbiowrapper", "cbiowrap");
        parentJob.addFile(cbioWrapDirGzSqw);
    }
    
    private Job runCbioWrap(String combinedMaf, String segFile, String gepFile, String cbiowrapDir) {
        if (gepFile == null){
            gepFile = "blank";
        }
        if (segFile == null){
            segFile = "blank";
        }
        if (combinedMaf == null){
            combinedMaf = "blank";
        }
        Job cbioWrap = getWorkflow().createBashJob("cbioWrap");
        Command cmd = cbioWrap.getCommand();
        cmd.addArgument(this.rPath + "/bin/Rscript " + this.cbioWrapper);
        cmd.addArgument(this.cbioWrapPath);
        cmd.addArgument(cbiowrapDir);
        cmd.addArgument(combinedMaf);
        cmd.addArgument(segFile);
        cmd.addArgument(gepFile);
        cmd.addArgument(this.hotspotGenesFile);
        cmd.addArgument(this.oncoKBFile);
        cmd.addArgument(this.ensembleConversonFile);
        cmd.addArgument(this.blackList);
        // set additional 
        cmd.addArgument("tar -zcvf " + cbiowrapDir + ".tar.gz" + " " + cbiowrapDir);
        cbioWrap.setMaxMemory(Integer.toString(this.cbiowrapMem * 1024));
        cbioWrap.setQueue(getOptionalProperty("queue", ""));
        return cbioWrap;
    }  
    
    private Job linkMAFs(String mafDir, String inMAFs, String combinedMaf) {
        Job linkMAF = getWorkflow().createBashJob("makeMafDir");
        String[] mafFiles = inMAFs.split(",");
        Command cmd = linkMAF.getCommand();
        cmd.addArgument("mkdir -p " + mafDir + ";\n");
        for (String mafFile : mafFiles){
            cmd.addArgument("ln -s " + mafFile + " " + mafDir);
        }
        cmd.addArgument("cd " + mafDir);
        cmd.addArgument("zcat $(ls *.maf.gz | head -1) | awk 'NR == 2' > " + combinedMaf);
        cmd.addArgument("for i in $(ls *.maf.gz); do zcat $i | awk 'NR > 2' >> " + combinedMaf + ";done");
        // set additional 
        linkMAF.setMaxMemory(Integer.toString(this.cbiowrapMem * 1024));
        linkMAF.setQueue(getOptionalProperty("queue", ""));
        return linkMAF;
    } 
    
    private Job copySegProcess(String inCopySegs, String combinedSeg) {
        Job combineCopySeg = getWorkflow().createBashJob("combineSegFiles");
        String[] copySegFiles = inCopySegs.split(",");
        Command cmd = combineCopySeg.getCommand();
        for (String segFile : copySegFiles){
            String tmpFL = this.tmpDir + getSampleName(segFile, ".seg") + ".txt";
            cmd.addArgument("head -1 " + segFile + " > " + tmpFL);
            if (segFile.contains(".sorted.filter.deduped.realign.recal")){
                cmd.addArgument("awk -F \"\\t\" -v OFS=\"\\t\" -v j=" + 
                    getSampleName(segFile, ".sorted.filter.deduped.realign.recal.seg") + 
                    " '{ if (NR>1) {print j,\"chr\"$2,$3,$4,$5,$6} }' " + segFile +
                    " >> " + tmpFL);
            } else{
                cmd.addArgument("awk -F \"\\t\" -v OFS=\"\\t\" -v j=" + 
                    getSampleName(segFile, ".seg") + 
                    " '{ if (NR>1) {print j,$2,$3,$4,$5,$6} }' " + segFile +
                    " >> " + tmpFL);
            }
        }
        cmd.addArgument("echo -e \"ID\\tchrom\\tloc.start\\tloc.end\\tnum.mark\\tseg.mean\" > " + combinedSeg + "\n"); 
        cmd.addArgument("cat $TEMP/*.txt | grep -v \"chrom\" >> " + combinedSeg);
        // set additional 
        combineCopySeg.setMaxMemory(Integer.toString(this.cbiowrapMem * 1024));
        combineCopySeg.setQueue(getOptionalProperty("queue", ""));
        return combineCopySeg;
    } 
    
    
    private Job postProcessRSEM(String inRSEMs, String inSTARs, String postProcessedRSEM) {
        Job postProcessRSEMGeneCounts = getWorkflow().createBashJob("post_process_RSEM");
        Command cmd = postProcessRSEMGeneCounts.getCommand();
        Map<String, List<String>> map = this.getRsemStarMap(inRSEMs, inSTARs);
        for (String key: map.keySet()){
            String geneCount = this.tmpDir + key + ".count";
            String geneRcount = this.tmpDir + key + ".rcount";
            String geneTPM = this.tmpDir + key + ".tpm";
            String geneFPKM = this.tmpDir + key + ".fpkm";
           
            String geneFPKMs = getFiles().get("RSEM_" + key + "_FPKM").getProvisionedPath();
            String geneTPMs = getFiles().get("RSEM_" + key + "_TPM").getProvisionedPath();
            cmd.addArgument("echo \"" + key + "\"" + ">" + geneTPM + "; cut -f6 " + geneTPMs + " | awk 'NR>1' >> " + geneTPM + ";\n");
            cmd.addArgument("echo \"" + key + "\"" + ">" + geneFPKM + "; cut -f6 " + geneFPKMs + " | awk 'NR>1' >> " + geneFPKM + ";\n");
            // rcounts
            String geneCounts = getFiles().get("RSEM_" + key + "_COUNTS").getProvisionedPath();
            String rtab = getFiles().get("STAR_"+key).getProvisionedPath();
            cmd.addArgument("echo \"" + key + "\" > " + geneCount + ";");
            cmd.addArgument("cut -f5 " + geneCounts + " | awk 'NR>1' >> " + geneCount + ";");
            cmd.addArgument("echo \"" + key + "\" > " 
                    + geneRcount 
                    + ";");
            cmd.addArgument("awk 'NR>4 {if ($4 >= $3) print $4; else print $3}' " 
                    + rtab + " >> " + geneRcount + ";");
            cmd.addArgument("cp " + rtab + " " + this.tmpDir + ";");
        }
        cmd.addArgument("STARG=`ls " + this.tmpDir + "*.tab | head -1`;");
        cmd.addArgument("if [ ! -z $STARG ]; then awk 'NR>3 {print $1}' $STARG | sed \"s/N\\_ambiguous/gene\\_id/\" > " + this.tmpDir + "sgene; fi;");
        cmd.addArgument("paste " + this.tmpDir + "sgene " + this.tmpDir + "*.rcount > " + postProcessedRSEM + ";\n");
        cmd.addArgument("paste " + this.tmpDir + "genes " + this.tmpDir + "*.tpm > " + postProcessedRSEM.replace("_RCOUNT.txt", "_GENES_TPM.txt"));
        cmd.addArgument("paste " + this.tmpDir + "genes " + this.tmpDir + "*.fpkm > " + postProcessedRSEM.replace("_RCOUNT.txt", "_GENES_FPKM.txt"));
        cmd.addArgument("paste " + this.tmpDir + "genes " + this.tmpDir + "*.count > " + postProcessedRSEM.replace("_RCOUNT.txt", "_GENES_COUNT.txt"));
        postProcessRSEMGeneCounts.setMaxMemory(Integer.toString(this.cbiowrapMem * 1024));
        postProcessRSEMGeneCounts.setQueue(getOptionalProperty("queue", ""));
        return postProcessRSEMGeneCounts; 
    }
    
    private String getSampleName(String fileBaseName, String extn){
        String[] sampleBaseName = fileBaseName.split(extn);
        List<String> sampleTokens = new ArrayList<String>(Arrays.asList(sampleBaseName[0].split("_")));
        List<String> sampleDesc = new ArrayList<String>();
        int st = 2;
        if (extn.contains(".seg")){
            st = 0;
        }
        for (int i = st; i < sampleTokens.size(); i++){
            String token = sampleTokens.get(i);
            sampleDesc.add(token);
        }
        String sampleName = String.join("_", sampleDesc);
        return sampleName;
    }
    
    private String getSTARMap(String rsemSampleName, String commaSeparatedSTAR){
        String starMap = new String (); 
        String[] starFilePaths = commaSeparatedSTAR.split(",");
        for (String starFile : starFilePaths){
            String starBaseName = FilenameUtils.getBaseName(starFile);
            String starSampleName = this.getSampleName(starBaseName, ".ReadsPerGene.out");
            if (starSampleName.equals(rsemSampleName)){
                starMap = starFile;
            }
        }
        return starMap;
    }
    
    private Map<String, List<String>> getRsemStarMap(String commaSeparatedRSEM, String commaSeparatedSTAR){
        /**
         * Given a list of comma Separated RSEM files;
         * get matching STAR
         */
        String[] rsemFilePaths = commaSeparatedRSEM.split(",");
        
        Map<String,List<String>> rsemStarMap = new HashMap<String, List<String>>();
        for (String rsemFile : rsemFilePaths){
            String rsemBaseName = FilenameUtils.getBaseName(rsemFile);
            String rsemSampleName = this.getSampleName(rsemBaseName, ".genes");
            List<String> vls = new ArrayList<String> ();
            vls.add(rsemFile);
            String starMap = getSTARMap(rsemSampleName, commaSeparatedSTAR);
            if ((starMap.equals("")) || (starMap == null)){
                Log.info("Missing star input");
                continue;
            }
            vls.add(starMap);
            if (vls.size() != 2 ){
                Log.debug("Skipping "+rsemSampleName+ "Missing STAR");
                continue;
            }
            rsemStarMap.put(rsemSampleName, vls);
        }
        return rsemStarMap;
    }
    
}