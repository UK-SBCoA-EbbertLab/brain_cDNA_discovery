// Import Modules
include {BASECALL; GATHER_BASECALL} from '../modules/basecall'

workflow BASECALLER {
        
    take:
        fast5_dir
        config
        id

    main:


    BASECALL(fast5_dir, config)
    GATHER_BASECALL(id, BASECALL.out.fastq.collect(), BASECALL.out.txt.collect())

}
