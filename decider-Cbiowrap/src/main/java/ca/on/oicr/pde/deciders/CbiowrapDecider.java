package ca.on.oicr.pde.deciders;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author prath@oicr.on.ca
 */
public class PureCNDecider extends OicrDecider {

    private final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;

    private String templateType = "EX";
    private String queue = "";
   
    private String bias="/.mounts/labs/PDE/Modules/sw/cnvkit/0.9.3/reference/mapping_bias_agilent_v6_hg19.rds";
    private String simp="/.mounts/labs/PDE/Modules/sw/cnvkit/0.9.3/reference/hg19_simpleRepeats.bed";
    private String hgbuild="hg19";
    private String pureCNMem= "32";
    
    
    
    
    private HashMap <String, String> fileMap = new HashMap<String, String>();

    private final static String[] VCF_GZ_METATYPES = {"application/vcf-gz", "application/vcf-4-gzip"};
    private final static String[] CALL_STATS_METATYPES = {"text/plain", "application/txt-gz"};
    private final static String MODEL_FIT_METATYPE = "application/tar-gzip";
    private final static String[] METATYPES = (String[])ArrayUtils.addAll(VCF_GZ_METATYPES, CALL_STATS_METATYPES, MODEL_FIT_METATYPE);
    private String[] allowedExtensions = new String[]{".tar.gz", ".muTect.tumor_only.snvs.vcf.gz", ".muTect.tumor_only.snvs.out", ".muTect.tumor_only.snvs.out.gz"};
    

    public PureCNDecider() {
        super();
        fileSwaToSmall = new HashMap<String, BeSmall>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("template-type", "Required. Set the template type to limit the workflow run "
                + "so that it runs on data only of this template type").withRequiredArg();
        parser.accepts("queue", "Optional: Set the queue (Default: not set)").withOptionalArg();
        parser.accepts("bias", "Optional:Provide the path for mapping bias rds. Default " + this.bias).withOptionalArg();
        parser.accepts("simp", "Optional:Provide the path for simple repeats bed file. Default " + this.simp).withOptionalArg();
        parser.accepts("hgbuild", "Optional:Provide the hgbuild version. Default " + this.hgbuild).withOptionalArg();
        parser.accepts("purecn-mem", "Optional:Provide memory requirements for pureCN workfflow. Default " + this.pureCNMem).withOptionalArg();
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setMetaType(Arrays.asList(METATYPES));
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
            this.templateType = options.valueOf("template-type").toString();
            if (!this.templateType.equals("EX")) {
                Log.warn("This workflow run may not schedule for --template-type " + this.templateType);
            }
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
        for (String p : filePaths) {
            if (p.endsWith("vcf.gz")){
                fileMap.put("VCF", p);
            }
            if (p.endsWith(".out") || p.endsWith(".out.gz")){
                fileMap.put("CALL", p);
            } 
            if (p.endsWith("tar.gz")){
                fileMap.put("MODEL-FIT", p);
            }
        }
        return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);
        String currentTtype = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        String currentTissueType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_tissue_type");
        String currentWorkflowName = returnValue.getAttribute(Header.WORKFLOW_NAME.getTitle());
        String fileMetaType = fm.getMetaType();

        if (null == currentTissueType) {
            return false; // we need only those which have their tissue type set
        }
        
        if (!Arrays.asList(METATYPES).contains(fileMetaType)){
           return false;
        } else {
            if (fileMetaType.equals(MODEL_FIT_METATYPE)){
                if (!currentWorkflowName.equals("CNVkit")){
                    return false;
                }
            }
        }

        // Filter the data of a different template type if filter is specified
        if (!this.templateType.equalsIgnoreCase(currentTtype)) {
            Log.warn("Excluding file with SWID = [" + returnValue.getAttribute(Header.FILE_SWA.getTitle())
                    + "] due to template type/geo_library_source_template_type = [" + currentTtype + "]");
            return false;
        }


        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        Log.debug("Number of files from file provenance = " + vals.size());

        // get files from study
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        // Override the supplied group-by value
        for (ReturnValue currentRV : vals) {
            boolean metatypeOK = false;
            boolean fileExtnOK = false;

            for (int f = 0; f < currentRV.getFiles().size(); f++) {
                try {
                    if (Arrays.asList(METATYPES).contains(currentRV.getFiles().get(f).getMetaType())) {
                        metatypeOK = true;
                    }
                    
                    fileExtnOK = this.identifyFilePath(currentRV.getFiles().get(f).getFilePath());
                    
                } catch (Exception e) {
                    Log.stderr("Error checking a file");
                }
            }

            if (!metatypeOK) {
                continue; // Go to the next value
            } 
            
            if (!fileExtnOK) {
                continue;
            }

            BeSmall currentSmall = new BeSmall(currentRV);
            fileSwaToSmall.put(currentRV.getAttribute(groupBy), currentSmall);

            String fileDeets = currentSmall.getIusDetails();
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
            String currVal = fileSwaToSmall.get(r.getAttribute(Header.FILE_SWA.getTitle())).getGroupByAttribute();
            List<ReturnValue> vs = map.get(currVal);
            if (vs == null) {
                vs = new ArrayList<ReturnValue>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }

        return map;
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
        Map<String, String> iniFileMap = super.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        String inVCF = this.fileMap.get("VCF");
        String outputFileNamePrefix = getExternalName(inVCF);
        String inCallStats = this.fileMap.get("CALL");
        String inModelFit = this.fileMap.get("MODEL-FIT");
        iniFileMap.put("input_vcf", inVCF);
        iniFileMap.put("call_stats_file", inCallStats);
        iniFileMap.put("cnvkit_model_fit", inModelFit);
        iniFileMap.put("output_filename_prefix", outputFileNamePrefix);
        iniFileMap.put("bias", this.bias);
        iniFileMap.put("simp", this.simp);
        iniFileMap.put("hgbuild", this.hgbuild);
        iniFileMap.put("purecn_mem", this.pureCNMem);
        return iniFileMap;
                

//        return super.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);
    }

    public static void main(String args[]) {

        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(PureCNDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }
    
    private boolean identifyFilePath(String filePath) {
        return Arrays.stream(allowedExtensions).anyMatch(entry -> filePath.endsWith(entry));
    }

    private String getExternalName(String inVCF) {
        String baseName = FilenameUtils.getBaseName(inVCF);
        String[] bNames = baseName.split("\\.");
        String sampleName = bNames[0];
        return sampleName;
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
        private String workflowName = null;

        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            workflowName = rv.getAttribute(Header.WORKFLOW_NAME.getTitle());
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode();
            tissueType = fa.getLimsValue(Lims.TISSUE_TYPE);
            extName = rv.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_external_name");
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
            groupByAttribute = iusDetails + ":" + groupID + ":" + extName;
            path = rv.getFiles().get(0).getFilePath() + "";
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
            return iusDetails;
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
    }

    /**
     * Utility function
     *
     * @param path
     * @param extension
     *
     * @return
     */
    public static String makeBasename(String path, String extension) {
        return path.substring(path.lastIndexOf("/") + 1, path.lastIndexOf(extension));
    }
}
