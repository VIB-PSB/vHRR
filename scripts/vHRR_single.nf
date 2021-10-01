#!/usr/bin/env nextflow

/*
 * A nextflow pipeline for the co-expression vHRR AFP protocol.
 */


// Scripts
script_prepare_starter_pack = "${params.scripts}/vHRR_prepare.py"
script_scan_HRR = "${params.scripts}/vHRR_scan.py"
script_final_run = "${params.scripts}/vHRR_collect_single.py"


/*
 * Step 1. Prepare starter-pack 
 * ----------------------------
 */

process prepare_starter_pack {

    module 'python'

    output:
    file 'go_terms.li' into go_terms_list

    script:
    """
        mkdir -p ${params.wkdir}/tmp
        mkdir -p ${params.wkdir}/tmp/enrichments
        mkdir -p ${params.out_folder}

        OMP_NUM_THREADS=1 python3 ${script_prepare_starter_pack} ${params.atlas} \
                                                                 ${params.go_gene} \
                                                                 ${params.wkdir}/tmp/starter_pack.npz \
                                                                 ${params.sample_weighing} \
                                                                 > go_terms.li
    """
}


go_terms_list.splitText() { it.trim() }
             .set{ go_terms_list_split }



/*
 * Step 2. Scan HRR range per GO 
 * -----------------------------
 */

process scan_HRR {

    module 'python'

    input:
    val go_term from go_terms_list_split

    output:
    file '.scan_dummy' into scan_dummy

    script:
    """
        touch .scan_dummy
        OMP_NUM_THREADS=1 python3 ${script_scan_HRR} ${params.atlas} \
                                                     ${params.wkdir}/tmp/starter_pack.npz \
                                                     ${go_term} \
                                                     ${params.sample_weighing} \
                                                     -out_folder ${params.wkdir}/tmp/enrichments \
                                                     > ${params.wkdir}/tmp/${go_term}_hrr.tmp
    """
}



/* Step 3. evaluate predictions 
 * -----------------------
 */

process evaluate_predictions {

    module 'python'

    input:
    file '.scan_dummy' from scan_dummy.collect()

    script:
    """
        cat ${params.wkdir}/tmp/*_hrr.tmp > ${params.wkdir}/tmp/hrr_sizes.tsv

        OMP_NUM_THREADS=1 python3 ${script_final_run} ${params.atlas} \
                                                      ${params.wkdir}/tmp/starter_pack.npz \
                                                      -hrr ${params.wkdir}/tmp/hrr_sizes.tsv \
                                                      -enrichments ${params.wkdir}/tmp/enrichments \
                                                      -all_features \
                                                      > ${params.out_folder}/predictions_single.tsv
    """
}
