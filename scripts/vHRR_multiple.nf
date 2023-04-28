#!/usr/bin/env nextflow

/*
 * A nextflow pipeline that fully automates itegrated co-expression gene function prediction.
 */


// Scripts
script_prepare_starter_pack = "${params.scripts}/vHRR_prepare.py"
script_scan_HRR = "${params.scripts}/vHRR_scan.py"
script_final_run = "${params.scripts}/vHRR_collect_multiple.py"



// Create channel for input atlasses
expression_compendia = Channel.fromPath(params.atlases)
expression_compendia.into { expression_compendia_1; 
                            expression_compendia_2 
                        }


/*
 * Step 1. Prepare starter-packs 
 * -----------------------------
 */

process prepare_starter_pack {

    module 'python'

    input:
    file atlas from expression_compendia_1

    output:
    file 'go_terms.li' into go_terms_list

    script:
    """
        mkdir -p ${params.tmpdir}
        mkdir -p ${params.tmpdir}/${atlas.baseName}
        mkdir -p ${params.tmpdir}/${atlas.baseName}/enrichments

        OMP_NUM_THREADS=1 python3 ${script_prepare_starter_pack} ${atlas} \
                                                                 ${params.go_gene} \
                                                                 ${params.tmpdir}/${atlas.baseName}/starter_pack.npz \
                                                                 > go_terms.li
    """
}


go_terms_list.last().splitText() { it.trim() }
                    .set{ go_terms_list_split }

atlases_terms = expression_compendia_2.combine(go_terms_list_split)



/*
 * Step 2. Scan HRR range per GO per atlas 
 * ---------------------------------------
 */

process scan_HRR {

    module 'python'

    input:
    tuple file(atlas), val(go_term) from atlases_terms

    output:
    file '.scan_dummy' into scan_dummy

    script:
    """
        touch .scan_dummy
        OMP_NUM_THREADS=1 python3 ${script_scan_HRR} ${atlas} \
                                                     ${params.tmpdir}/${atlas.baseName}/starter_pack.npz \
                                                     ${go_term} \
                                                     -out_folder ${params.tmpdir}/${atlas.baseName}/enrichments \
                                                     > ${params.tmpdir}/${atlas.baseName}/${go_term}_hrr.tmp
    """
}



/*
 * Step 3. Assigning GO terms
 * --------------------------
 */

process evaluate_predictions {

    module 'python'

    input:
    file '.scan_dummy' from scan_dummy.collect()

    script:
    """
        cat ${params.tmpdir}/*/*_hrr.tmp > ${params.tmpdir}/hrr_sizes.tsv

        OMP_NUM_THREADS=1 python3 ${script_final_run} ${params.go_gene} \
                                                      ${params.tmpdir} \
                                                      -hrr ${params.tmpdir}/hrr_sizes.tsv \
                                                      > ${params.out_folder}/predictions_multiple.tsv
    """
}
