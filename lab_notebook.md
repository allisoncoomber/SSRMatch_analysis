# Evaluating multiclass classifier for P. infestans

September 9, 2020

I want to test my classifier for P. infestans lineages to get some measures of how well it is working. I am planning on using several binary classifier evaluation metrics for each class individually. This website follows a similar protocol and has relevant explanations [helpful site](https://parasite.id/blog/2018-12-13-model-evaluation/). I have read about these methods elsewhere in the literature but I find this site to be particularly clear. 

So the dataset for my classifier consists of samples which have all 12 loci present and samples which are missing some loci. I have decided to use only the samples which have all 12 loci present for the classifier reference set for now. To test the classifier I am subsampling this reference set to create five unqiue samples and corresponding sample reference sets (where the corresponding sample reference set is the entire reference set minus the sample). I wrote an R script to do this five times making sure the samples do not overlap. [Script here](https://github.com/allisoncoomber/SSRMatch_analysis/blob/master/generate_samples.R)

I then launch the app and run a sample through it as the test data with the "reference set" as the sample reference set corresponding to the sample. This has to be repeated five times, one for each sample. I am having difficulty with R picking up the delimiters in a CSV file (I think because the Source column contains commas itself) so I convert the output from my above script to tab delimited text before feeding it into the app. The applet code to run a sample against it's reference set for example is:
```
library(shinytest)
library(shiny)
library (poppr)
library (reshape2)
library (dplyr)
library(magrittr)

BRUVOS_MATCH_THRESHOLD <-  0.09929323
K <- 1

# TODO: Low Priority: Consider setting some columns as factors
DATABASE <- read.table("/Users/Allison/Documents/NCSU/Ristaino Lab/Workshop/sample_5_reference.txt", header=TRUE, encoding="UTF-8",stringsAsFactors = FALSE, colClasses = c(rep("character",20)))
DESCRIPTIONS <- read.table("descriptions.txt", header=FALSE, sep="\t", stringsAsFactors = FALSE)

 load_user_csv <- function(userCSV) {
   tryCatch(read.csv(userCSV, 
                     strip.white = TRUE, stringsAsFactors = FALSE, colClasses = c(rep("character",20)),
                     header=TRUE), warning = function(cond) {
                       output_warn <- print("We were unable to parse your CSV file. \n
                                            This could be caused by a few different things:
                                            \t - Check for extraneous commas in your file.
                                            \t - Make sure there is a new line at the end of your file. 
                                            \t   To do this, open your CSV in a text editor and hit enter 
                                            \t   (return) after the last line, then save your file and 
                                            \t   try again.")
                       validate(
                         need(output_warn == "", print(output_warn))
                       )
                     }
   )
 }
  
 format_input_for_bruvos_calculation <- function(data) {
   metadata <- select(data, "Sample","Genotype","Region","Country","Date","Host","Missing", "Present")
   tryCatch(
    data %>%
      select("D13","PinfSSR8","PinfSSR4","Pi04","Pi70","PinfSSR6","Pi63","PiG11","Pi02","PinfSSR11","PinfSSR2","Pi4B") %>%
      df2genind(ploidy=3, sep="/", ind.names=metadata[, 1], pop=metadata[, 2]), 
    warning = function(cond) {
     output_warn <- print("Ploidy is incorrect.")
     validate(
       need(output_warn == "", print(output_warn))
     )
    }
   )
 }

# A matrix of bruvos distances
get_bruvos <- function(query, reference) {
  withProgress(message = 'Calculating Bruvos Distance', value=0.4, detail = "This usually takes about 5 minutes, please be patient.", {

    return(
      bruvo.between(query, reference, replen=c(2,2,2,2,3,2,3,2,2,2,2,2), add=TRUE, loss=TRUE) %>%
      as.matrix() %>%
      melt(varnames = c("InputSample", "ReferenceSample"))
    )
  })
}

get_best_match <- function(bruvo_df, user_input) {
  withProgress(message = 'Calculating Bruvos Distance', value=1, detail = "This usually takes about 5 minutes, please be patient.", {
    return(
     bruvo_df %>%
        filter(InputSample %in% user_input$Sample) %>%
        filter(!(ReferenceSample %in% user_input$Sample)) %>%
        group_by(InputSample) %>%
        slice_min(value, n=K, with_ties=FALSE)
    )
  })
}

# TODO: Split into functions that obey single responsibility principle
clean_up_best_match_table <- function(bestmatch) {
  bestmatch %>%
    left_join(select(DATABASE, "Sample","Genotype","Region","Country","Date","Host","Missing", "Present"), by = c("ReferenceSample" = "Sample")) %>%
    rename(c("Input Sample"="InputSample", "Matching Genotype"="Genotype", "Bruvo's Distance"="value")) %>%
    select("Input Sample", "ReferenceSample", "Matching Genotype", "Region", "Country", "Date", "Host", "Bruvo's Distance")
}

weak_bruvos_matches <- function(bruvosDistances, threshold) {
  return(
    bruvosDistances %>%
      filter(`Bruvo's Distance` > threshold)
  )
}

strong_bruvos_matches <- function(bruvosDistances, threshold) {
  return(
    bruvosDistances %>%
      filter(`Bruvo's Distance` < threshold)
  )
}

as_no_matches_string <- function(bruvosDistances) {
  bruvos_strings <-
  bruvosDistances %>%
    mutate(no_match_string = paste("No close matches found for sample", `Input Sample`, ".")) # Paste the name and the rest of the string
  return(paste(bruvos_strings$no_match_string, sep = "\n", collapse = "\n")) # join the elements of the column with a newline
}



shinyServer(
  function(input, output) {
    calculation <- reactive({ 
      
      inFile <- input$SSRinput
      
      validate(
        need(!is.null(inFile), "Please select a data set.")
      )
        
      userInput <- load_user_csv(inFile$datapath)
      validate(
        need(identical(c("Sample","Genotype","Region","Country","Date","Host","D13","PinfSSR8","PinfSSR4","Pi04","Pi70","PinfSSR6","Pi63","PiG11","Pi02","PinfSSR11","PinfSSR2","Pi4B","Missing", "Present"), colnames(userInput)), "Header is incorrect.")
      )
      
      if (is.null(inFile))
        return(NULL)
        
      # Run calculation
      REFERENCE <- format_input_for_bruvos_calculation(DATABASE)
      validate(
          need(format_input_for_bruvos_calculation(userInput), "Ploidy is incorrect.")
      )
      query <- format_input_for_bruvos_calculation(userInput)
      return( 
        get_bruvos(query, REFERENCE) %>%
        get_best_match(userInput) %>%
        clean_up_best_match_table() %>%
          right_join(userInput, by = c("Input Sample" = "Sample")) %>%
          select("Input Sample", "Matching Genotype", "Genotype", "Bruvo's Distance")
      )
    })
    
    output$closestmatches <- renderTable({calculation() %>%
        strong_bruvos_matches(BRUVOS_MATCH_THRESHOLD)})
    
    output$no_match <- renderText({calculation() %>% 
        weak_bruvos_matches(BRUVOS_MATCH_THRESHOLD) %>% 
        as_no_matches_string}) 

    
    output$detail <- renderTable({
      # Table that has columns: sample name, description
      # The sample name comes from: calculation(), a reactive variable, in that object the column
      #   is called `Input Sample`. `Matching_Genotype`
      # description comes from descriptions df, matching column is `V1` in descriptions.
      # We are matching genotype from these two tables. 
      
      if (is.null(calculation())) {
        return(NULL)
      } else {
        return(
          DESCRIPTIONS %>%
            # Combine columns
            mutate(Description = paste(V3, V4, sep=' ')) %>%
            # Join tables based on matching genotypes
            inner_join(calculation() %>% strong_bruvos_matches(BRUVOS_MATCH_THRESHOLD),
                       by = c("V1"="Matching Genotype"), copy=TRUE) %>%
            # Subset table to get only rows of interest
            select(`Input Sample`, `Description`)
        )
      }
    })
  })
    
```


When I put a sample into the app to run it I first need to remove the numerical ID columns at the beginning and end. I was originally using these to make sure I was getting unique random samples as desired in the previous step. I then copy and pasted the output table of matches into an tab delimited text file and named it "K1_sampleX_results.txt" where X is the sample number. The K1 here represents that I am looking only at the one closest match. 

Next, I wrote a little script to make a "true/false" determination for the genotypes as well as generate a confusion matrix for each sample run. [Script here](https://github.com/allisoncoomber/SSRMatch_analysis/blob/master/TF_matches.R)

Going forward I plan to sum these confusion matrices all into one confusion matrix and then generate a bunch of summary statistics based on the ratio of True Positives, False Positives, True Negatives, and False Negatives for each class. 

***
September 10, 2020

Today I summarized the confusion matrices I generated yesterday into one large matrix using R. [Script here](https://github.com/allisoncoomber/SSRMatch_analysis/blob/master/summary_statistics_SSR_matcher.R)

I then reformatted this table in Excel and added in 0s for any lineages which did not have any matches. Actual matches are across the top and predicted matches are across the y-axis. 


![Confusion Matrix](https://github.com/allisoncoomber/SSRMatch_analysis/blob/master/Images/confusion_matrix.png)
