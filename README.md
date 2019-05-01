# somatic-variant-caller

An example run of the pipeline may look something like this:
```bash
nextflow run lifebit-ai/somatic-variant-caller --genome hg19 --samples s3://lifebit-featured-datasets/pipelines/gatk-somatic-data/samples.tsv --bed s3://lifebit-featured-datasets/pipelines/DeepVariant/quick_test/test_nist.b37_chr20_100kbp_at_10mb.bed
```

## Dependencies 
[Nextflow](https://www.nextflow.io/)
[Docker](https://www.docker.com/)

## Deploit
The pipeline can also be run on [Deploit](https://lifebit.ai/deploit), a bioinformatics platform to help you run your analysis over the cloud which is free for indiviudal users

![deploit](https://raw.githubusercontent.com/lifebit-ai/ecw-converter/master/images/deploit.png)

See an example job run on the Deploit platform [here](https://deploit.lifebit.ai/public/jobs/5cc97409b1d9e100b2672b33)
