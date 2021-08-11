# rcea
This is the repository for the `rcea` package, which accompanies a short course on model-based cost-effectiveness analysis (CEA) with `R`. A range of models are covered including time-homogeneous and time-inhomogeneous Markov cohort models, partitioned survival models, and semi-Markov individual patient simulations. In addition, the course shows how simulated costs and QALYs from a probabilistic sensitivity analysis can be used for decision-analysis within a cost-effectiveness framework. Analyses are conducted using both base `R` and the `R` package [hesim](https://hesim-dev.github.io/hesim/).

The course materials are available at https://hesim-dev.github.io/rcea.

## Installation and setup
All required `R` packages and course materials can be installed with the following steps.

1. Open an R session. We recommend using [RStudio](https://rstudio.com/).

2. Install the `rcea` package from GitHub, which will also install all other required packages.

    ```r
    # install.packages("devtools") # You must install the "devtools" R package first.
    devtools::install_github("hesim-dev/rcea")
    ```

3. Create a new project in your desired directory. 

    ```r
    # Create a project named "rcea-exercises" within a directory named "Projects"
    usethis::create_project("~/Projects/rcea-exercises") 
    ```

4. Add the course materials (`R` scripts for the tutorials) to your new project.

    ```r
    rcea::use_rcea("~/Projects/rcea-exercises")
    ```

## Tutorials
The course contains six tutorials:

1. **Markov Cohort Model**: A simple time-homogeneous Markov cohort model with fixed parameter values.

2. **Incorporating Probabilistic Sensitivity Analysis (PSA)**: The Markov cohort model is re-analyzed using suitable probability distributions for the parameters. 

3. **Markov Cohort Model with hesim**: The second tutorial---programmed primarily using base `R`--- is repeated using the `R` package `hesim`. 

4. **Semi-Markov Multi-state state Model**: A semi-Markov multi-state model is fit to patient-level data and outcomes are simulated using an individual patient simulation. 

5. **Partitioned Survival Model**: The data from the fourth tutorial is refit using partitioned survival analysis and state probabilities are computed using the "area under the curve" technique.

6. **Cost-effectiveness Analysis (CEA)**: CEA is performed using the cost and QALY output of the PSA from the fourth tutorial. A number of methods are used to represent decision uncertainty (e.g. cost-effectiveness planes, cost-effectiveness acceptability curves, and cost-effectiveness acceptability frontiers), and value of information analysis is conducted. 

## Learning R
For those new to `R`, we recommend the following free online resources:

* [R for Data Science](https://r4ds.had.co.nz/) teaches `R` for data science with the [`tidyverse`](https://www.tidyverse.org/).

* [An introduction to `R`](https://cran.r-project.org/doc/manuals/r-devel/R-intro.pdf) is official [CRAN](https://cran.r-project.org/) documentation covering foundational concepts and use of base `R`. 

A list of additional resources is also available [here](https://stackoverflow.com/tags/r/info).

## data.table
We use [`data.table`](https://rdatatable.gitlab.io/data.table/) to summarize output because it is very fast when working with large datasets, as is often produced by simulation models. For those more familiar with [`dplyr`](https://dplyr.tidyverse.org/), a nice comparison between `dplyr` and `data.table` can be found [here](https://atrebas.github.io/post/2019-03-03-datatable-dplyr/).



