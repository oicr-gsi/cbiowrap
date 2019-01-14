package ca.on.oicr.pde.deciders;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;

/**
 *
 * @author prath@oicr.on.ca
 */
public class CbiowrapDecider extends OicrDecider {

    private final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;

    private String templateType = "EX,WT";
    private String[] allowedTemplateTypes = {"EX", "WT"};
    private String queue = "";
    private String[] allowedExtensionTypes = {".ReadsPerGene.out.tab", ".genes.results", ".maf.txt.gz", ".seg"};
    private String studyName;
    private String studyTitle;
    private String hotSpotGenes = "/.mounts/labs/TGL/gsi/databases/20161121_Chang_hotspot_list.tsv";
    private String ensFile = "/.mounts/labs/TGL/gsi/databases/ensemble_conversion.txt";
    private String oncoKB = "/.mounts/labs/TGL/gsi/databases/20170412_oncoKB.tsv";
    private String blackList = "/.mounts/labs/TGL/gsi/databases/blacklist.txt";

    private final static String TXT_METATYPE = "text/plain";
    private final static String TXT_GZ_METATYPE = "application/txt-gz";
    private final String[] allowedMetatypes = {"text/plain", "application/txt-gz"};
//    private String tumorType;
//    private List<String> results;
    //private String groupKey;

    public CbiowrapDecider() {
        super();
        fileSwaToSmall = new HashMap<String, BeSmall>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("template-type", "Required. Set the template type to limit the workflow run "
                + "so that it runs on data only of this template type").withRequiredArg();
        parser.accepts("queue", "Optional: Set the queue (Default: not set)").withRequiredArg();
        parser.accepts("study-name", "Required. Specify study name, e.g. TGL07, OCT").withRequiredArg();
        parser.accepts("template-type", "Optional. Specify template-type. Default WT,EX").withOptionalArg();
        parser.accepts("ensemble-file", "Optional. Specify ensemble gene text file").withOptionalArg();
        parser.accepts("hotspot-genes", "Optional. Specify ensemble gene text file").withOptionalArg();
        parser.accepts("blacklist-tsv", "Optional. Specify ensemble gene text file").withOptionalArg();
        parser.accepts("onco-kb", "Optional. Specify ensemble gene text file").withOptionalArg();
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setMetaType(Arrays.asList(allowedMetatypes));
        this.setHeadersToGroupBy(Arrays.asList(Header.FILE_SWA));

        ReturnValue rv = super.init();
        rv.setExitStatus(ReturnValue.SUCCESS);

        //Group by sample if no other grouping selected
        if (this.options.has("group-by")) {
            Log.error("group-by parameter passed, but this decider does not allow overriding the default grouping (by Donor + Library Type)");
        }

        if (this.options.has("queue")) {
            this.queue = options.valueOf("queue").toString();
        }

        if (this.options.has("template-type")) {
            String enteredTemplateTypes = options.valueOf("template-type").toString();
            String[] inputTT = enteredTemplateTypes.split(",");
            for (String tt : inputTT) {
                if (!Arrays.asList(this.allowedTemplateTypes).contains(tt)) {
                    Log.warn("Template type not in the list of allowed templates. Works best for template types WT and EX");
                }
            }
            this.templateType = enteredTemplateTypes;
            if (!this.templateType.equals("EX,WT")) {
                Log.warn("NOTE SPECIFY BOTH template-types SUPPORTED, WE CANNOT GUARANTEE MEANINGFUL RESULTS WITH OTHER TEMPLATE TYPES");
            }
        }
        if (this.options.has("study-name")) {
            if (!options.hasArgument("study-name")) {
                Log.error("--study-name requires study title, e.g. OCT, TGL07");
                rv.setExitStatus(ReturnValue.INVALIDARGUMENT);
                return rv;
            } else {
                this.studyTitle = options.valueOf("study-name").toString();
            }
        }
        
        if (this.options.has("ensemble-file")) {
            this.ensFile = options.valueOf("ensemble-file").toString();
        }
        if (this.options.has("onco-kb")) {
            this.oncoKB = options.valueOf("onco-kb").toString();
        }
        if (this.options.has("hotspot-genes")) {
            this.hotSpotGenes = options.valueOf("hotspot-genes").toString();
        }
        if (this.options.has("blacklist-tsv")) {
            this.blackList = options.valueOf("blacklist-tsv").toString();
        }
        
        return rv;
    }

    /**
     * Final check
     *
     * @param commaSeparatedFilePaths
     * @param commaSeparatedParentAccessions
     *
     * @return
     */
    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        String[] filePaths = commaSeparatedFilePaths.split(",");
        boolean haveFiles = false;
        List<String> studyNames = new ArrayList<String>();
        for (String p : filePaths) {
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p)) {
                    continue;
                }
                String tt = bs.getTissueType();

                if (!tt.isEmpty()) {
                    haveFiles = true;
                    String sn = bs.getStudyTitle();
                    studyNames.add(sn);
                }
            }
        }
        if (haveFiles) {
//            this.studyName = 
            HashSet<String> hsetStudyNames = new HashSet(studyNames);
            Iterator<String> itr = hsetStudyNames.iterator();
            this.studyName = itr.next().toString();
            return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        }
       Log.error("Data not available, WON'T RUN");
        return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);
        String currentTtype = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        String currentTissueType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_tissue_type");

        if (null == currentTissueType) {
            return false; // we need only those which have their tissue type set
        }

        // Filter the data of a different template type if filter is specified
        if (!Arrays.asList(this.allowedTemplateTypes).contains(currentTtype)) {
            Log.warn("Excluding file with SWID = [" + returnValue.getAttribute(Header.FILE_SWA.getTitle())
                    + "] due to template type/geo_library_source_template_type = [" + currentTtype + "]");
            return false;
        }


        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        // get files from study
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        // Override the supplied group-by value
        for (ReturnValue currentRV : vals) {
            boolean metatypeOK = false;
            boolean fileExtnOK = false;

            for (int f = 0; f < currentRV.getFiles().size(); f++) {
                try {
                    if (currentRV.getFiles().get(f).getMetaType().equals(TXT_METATYPE)) {
                        metatypeOK = true;
                    }
                    if (currentRV.getFiles().get(f).getFilePath().endsWith(".ReadsPerGene.out.tab") 
                            || currentRV.getFiles().get(f).getFilePath().endsWith(".genes.results")
                            || currentRV.getFiles().get(f).getFilePath().endsWith(".seg")
                            || currentRV.getFiles().get(f).getFilePath().endsWith(".maf.txt.gz")){
                        fileExtnOK = true;
                    } else {
//                        Log.debug("Undesired file type "+currentRV.getFiles(). get(f).getFilePath());
                        continue;
                    }
                } catch (Exception e) {
                    Log.stderr("Error checking a file");
                }
            }
            
            if (metatypeOK && fileExtnOK) {
                //
            } else {
                continue;
            }

            BeSmall currentSmall = new BeSmall(currentRV);
            fileSwaToSmall.put(currentRV.getAttribute(groupBy), currentSmall);
            //make sure you only have the most recent single file for each
            //sequencer run + lane + barcode + workflow-name
            String fileDeets = currentSmall.getIusDetails() + "_" + currentSmall.getWorkflowDetails();
            Date currentDate = currentSmall.getDate();

            //if there is no entry yet, add it
            if (iusDeetsToRV.get(fileDeets) == null) {
                iusDeetsToRV.put(fileDeets, currentRV);
            } //if there is an entry, compare the current value to the 'old' one in
            //the map. if the current date is newer than the 'old' date, replace
            //it in the map
            else {
                ReturnValue oldRV = iusDeetsToRV.get(fileDeets);
                BeSmall oldSmall = fileSwaToSmall.get(oldRV.getAttribute(Header.FILE_SWA.getTitle()));
                Date oldDate = oldSmall.getDate();
                if (currentDate.after(oldDate)) {
                    iusDeetsToRV.put(fileDeets, currentRV);
                }
            }
            
        }
        
        //only use those files that entered into the iusDeetsToRV
        //since it's a map, only the most recent values
        List<ReturnValue> newValues = new ArrayList<ReturnValue>(iusDeetsToRV.values());
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();

        //group files according to the designated header (e.g. sample SWID)
        for (ReturnValue r : newValues) {
//            String rootSampleName = fileSwaToSmall.get(r.getAttribute(Header.FILE_SWA.getTitle())).getRootSampleName();
            String currVal = fileSwaToSmall.get(r.getAttribute(Header.FILE_SWA.getTitle())).getGroupByAttribute();
            Log.debug(currVal);
            List<ReturnValue> vs = map.get(currVal);
            String filePath = r.getFiles().get(0).getFilePath();
            if (vs == null) {
                vs = new ArrayList<ReturnValue>();
            }
            
            if (filePath.endsWith(".genes.results") || 
                    filePath.endsWith(".ReadsPerGene.out.tab") ||
                    filePath.endsWith(".seg") ||
                    filePath.endsWith(".maf.txt.gz")){
                    vs.add(r); 
            } else {
//                Log.debug("Rejecting " + r.getFiles().get(0).getFilePath() );
                continue;
            }
            map.put(currVal, vs);
            //this.groupKey = currVal;
        }
        
        Map<String, List<ReturnValue>> filteredMap = new HashMap<String, List<ReturnValue>>();
        for (Entry<String, List<ReturnValue>> e :  map.entrySet()){
            String study = e.getKey();
            List<ReturnValue> mapValues = e.getValue();
       
            List<ReturnValue> finalRVs = new ArrayList<ReturnValue>();
            for (ReturnValue rV : mapValues) {
                String deets = fileSwaToSmall.get(rV.getAttribute(Header.FILE_SWA.getTitle())).getIusDetails();
                ArrayList<FileMetadata> filePaths = rV.getFiles();
                for (FileMetadata fm : filePaths) {
                    String fileName = fm.getFilePath();
                    if (fileName.endsWith(".genes.results")) {
                        Log.debug(deets + ":" + fileName);
                        ReturnValue rsemMapToDeets = rV;
                        ReturnValue starMapToDeets = getMappingSTARFileDeets(map, deets, study);
                        finalRVs.add(rsemMapToDeets);
                        finalRVs.add(starMapToDeets);
                    } else if (fileName.endsWith(".maf.txt.gz") || fileName.endsWith(".seg")){
                        finalRVs.add(rV);
                    }
                }
            }
            filteredMap.putIfAbsent(study, finalRVs);
        }
        return filteredMap;
    }
    
    protected ReturnValue getMappingSTARFileDeets(Map<String, List<ReturnValue>> map, String iusKey, String groupKey){
        ReturnValue starFilePath = new ReturnValue();
        List<ReturnValue> mapValues = new ArrayList<ReturnValue>(map.get(groupKey));
        for (ReturnValue rV : mapValues){
            String deets = fileSwaToSmall.get(rV.getAttribute(Header.FILE_SWA.getTitle())).getIusDetails();
            if (deets.equals(iusKey)){
                ArrayList<FileMetadata> filePaths = rV.getFiles();
                for (FileMetadata fm : filePaths){
                    String fileName = fm.getFilePath();
                    if (fileName.endsWith(".ReadsPerGene.out.tab")){
                        Log.debug(deets + ":" + fileName);
                        starFilePath = rV;
                    } else {
//                        Log.debug("Reading next ReturnValue rV for the .ReadsPerGene.out.tab file");
                        continue;
                    }
                }
            } else {
//                Log.debug("Skipping unmatched " + deets);
                continue;
            }
        }
        return starFilePath;
    }

    @Override
    protected String handleGroupByAttribute(String attribute) {
        String a = super.handleGroupByAttribute(attribute);
        BeSmall small = fileSwaToSmall.get(a);
        if (small != null) {
            return small.getGroupByAttribute();
        }
        return attribute;
    }

    

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {

        String[] filePaths = commaSeparatedFilePaths.split(",");
        List<String> rsemGeneCounts = new ArrayList<String>();
        List<String> starGeneCounts = new ArrayList<String>();
        List<String> mafFiles = new ArrayList<String>();
        List<String> segFiles = new ArrayList<String>();

        for (String p : filePaths) {

            for (BeSmall bs : fileSwaToSmall.values()) {
                String tt = bs.getTissueType();
                if (!tt.isEmpty()) {
                    if (!bs.getPath().equals(p)) {
                        continue;
                    }
                    if (p.endsWith(".genes.results")) {
                        rsemGeneCounts.add(p);
                    } else if (p.endsWith(".ReadsPerGene.out.tab")) {
                        starGeneCounts.add(p);
                    } else if (p.endsWith(".seg")) {
                        segFiles.add(p);
                    } else if (p.endsWith(".maf.txt.gz")) {
                        mafFiles.add(p);
                    }
                } else {
                    Log.error("THE DONOR does not have data to run the workflow");
                    abortSchedulingOfCurrentWorkflowRun();
                }
            }

        }
        Map<String, String> iniFileMap = super.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);

        iniFileMap.put("data_dir", "data");
        if (!rsemGeneCounts.isEmpty()){
            iniFileMap.put("input_rsem_counts", String.join(",", rsemGeneCounts));
        }
        if (!starGeneCounts.isEmpty()){
            iniFileMap.put("input_star_rtabs", String.join(",", starGeneCounts));
        }
        iniFileMap.put("input_mafs", String.join(",", mafFiles));
        iniFileMap.put("input_segs", String.join(",", segFiles));

        if (!this.queue.isEmpty()) {
            iniFileMap.put("queue", this.queue);
        }
        iniFileMap.put("study_title", this.studyName);
        iniFileMap.put("blacklist", this.blackList);
        iniFileMap.put("ensemble_gene_file", this.ensFile);
        iniFileMap.put("oncokb", this.oncoKB);
        iniFileMap.put("chang_hotspot", this.hotSpotGenes);

        //remove input_files, this is handled by rsem_inputs and star_inputs
        iniFileMap.remove("input_files");

        return iniFileMap;
    }

    public static void main(String args[]) {

        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(CbiowrapDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }
    private class BeSmall {

        private Date date = null;
        private String iusDetails = null;
        private String groupByAttribute = null;
        private String tissueType = null;
        private String path = null;
        private String extName = null;
        private String groupID = null;
        private String groupDescription = null;
        private String workflowDetails = null;
        private String sampleNameDetails = null;
        private String rootSampleName = null;
        private String studyTitle = null;

        
        
        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            studyTitle = rv.getAttribute(Header.STUDY_TITLE.getTitle());
            rootSampleName = rv.getAttribute(Header.ROOT_SAMPLE_NAME.getTitle());
            workflowDetails = rv.getAttribute(Header.WORKFLOW_NAME.getTitle());
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode(); //+ "_"+ workflowDetails;
            sampleNameDetails = fa.getLibrarySample() + "_" + fa.getSequencerRun() + fa.getLane() + fa.getBarcode();
            tissueType = fa.getLimsValue(Lims.TISSUE_TYPE);
            extName = rv.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_external_name");
            //fa.getLimsValue(Lims.TUBE_ID);
            if (null == extName || extName.isEmpty()) {
                extName = "NA";
            }
            groupID = fa.getLimsValue(Lims.GROUP_ID);
            if (null == groupID || groupID.isEmpty()) {
                groupID = "NA";
            }
            groupDescription = fa.getLimsValue(Lims.GROUP_DESC);
            if (null == groupDescription || groupDescription.isEmpty()) {
                groupDescription = "NA";
            }
            groupByAttribute = fa.getStudy() + ":" +  
                    fa.getMetatype() + ":" + workflowDetails;
            path = rv.getFiles().get(0).getFilePath() + "";
            
        }

        public String getRootSampleName() {
            return rootSampleName;
        }

        public String getSampleNameDetails() {
            return sampleNameDetails;
        }

        public String getWorkflowDetails() {
            return workflowDetails;
        }

        public Date getDate() {
            return date;
        }

        public void setDate(Date date) {
            this.date = date;
        }

        public String getGroupByAttribute() {
            return groupByAttribute;
        }

        public void setGroupByAttribute(String groupByAttribute) {
            this.groupByAttribute = groupByAttribute;
        }

        public String getTissueType() {
            return tissueType;
        }

        public String getIusDetails() {
            return this.iusDetails;
        }

        public void setIusDetails(String iusDetails) {
            this.iusDetails = iusDetails;
        }

        public String getPath() {
            return path;
        }

        public String getExtName() {
            return extName;
        }

        public String getGroupID() {
            return groupID;
        }

        public String getGroupDescription() {
            return groupDescription;
        }

        public void setPath(String path) {
            this.path = path;
        }
        public String getStudyTitle() {
            return studyTitle;
        }
    }


}
