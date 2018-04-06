# guidemapper-nf
Annotate CRISPR and shRNA screening libraries

## Installation
Install Nextflow following the instructions at https://www.nextflow.io/docs/latest/getstarted.html

The most convenient way is to install guidemapper-nf is to use nextflow's built-in 'pull' command
```bash
nextflow pull zuberlab/guidemapper-nf
```

## Test
Run the pipeline with test data
```bash
nextflow run zuberlab/guidemapper-nf --inputDir test_input --parameters test_parameters.txt
```

## Documentation
```bash
nextflow run zuberlab/guidemapper-nf --help
```

## Credits
Nextflow:  Paolo Di Tommaso - https://github.com/nextflow-io/nextflow
