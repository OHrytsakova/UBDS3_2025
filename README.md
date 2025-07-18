# UBDS3_2025
Teaching materials for the project "Community ecology based on environmental DNA data", which takes place July 28 - August 2, 2025 at the [Ukrainian Biological Data Science Summer School](https://www.bds3.org/) in Uzhhorod, Ukraine.

Environmental DNA offers a rapid, unbiased, and (moderately) low-cost approach for assessing species composition at a given site. However, interpreting large metabarcoding datasets might not be easy and requires not only substantial bioinformatics skills but also a deep understanding of the diversity and ecology of the focal group. In this project, we will dive into a real-world case of soil fungi data from Ukraine to explore both opportunities and pitfalls of such kind of data and find out what use we can make of it to answer the basic ecological questions.

## Repository structure
- **data**
    - ./lotus2_ITS2
        - ./ecm_physeq.Rdata
        - ./OTU.fna
    - ./lotus2_SSU_ASVs
        - ./amf_physeq.Rdata
        - ./OTU.fna
    - ./habitat_types.csv
    - ./tp-prylutskyi23.pdf
    - ./Ukraine_poly_adm0.Rdata
    - ./Ukraine_poly_adm1.Rdata
- **figures**
- **functions**
    - ./phyloseq_to_MDT_excel.R
- **scripts**
    - ./01_phyloseq.R
- UBDS3_2025.Rproj

Here you can see an example of how to deal with _phyloseq.R_ data object: 

[![Watch the video](https://img.youtube.com/vi/L8Orz34z_nk/0.jpg)](https://www.youtube.com/watch?v=L8Orz34z_nk)

### Useful links
1. [phyloseq: Explore microbiome profiles using R](https://joey711.github.io/phyloseq/index.html)

2. [Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses](https://f1000research.com/articles/5-1492/v2)

3. [Handling and analysis of high-throughput microbiome census data](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)

4. [Lab 7: Phyloseq](https://web.stanford.edu/class/bios221/labs/phyloseq/lab_phyloseq.html)

