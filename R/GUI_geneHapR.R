if(tools:::.OStype() == "windows"){
    build_app_geneHapR <- function(...) {
        shinyApp(ui = ui(), server = server)
    }
}

if(tools:::.OStype() == "windows"){
    startGUI.geneHapR <- function(...){
        app <- build_app_geneHapR()
        runApp(app)
    }
}

if(tools:::.OStype() == "windows"){
    p <-"
<b>TODO:</b>
 Reduce: ggpubr
 optimize: UI
 add: save as hmp
 add: filter largeVCF
 add: fasta2hap
 add: hmp2hap
 add: table2hap
 add: filter_hap
 add: adjust hapresult
 add: align and trim sequences
 add: vcf2hmp
 add: vcf2p.link
 add: site effect
 add: import mapdata base
<b>Done</b>
 Done: import fasta
 Done: import vcf
 Done: import p.link
 Done: filter vcf
 Done: filter p.link
 Done: vcf2hap
 Done: p.link2hap
 Done: seqs2hap
 Done: save hapResult to txt
 Done: save hapSummary to txt
 Done: visualize hap result as table
 Done: visualize variants on gene model
 Done: visualize haplonet
 Done: visualize geo-distribution
 Done: pheno comparison
"
    h_en <- '
<h2>Main steps:</h2>
<h3>1. Import data</h3>
<b>1.0 Basic steps:</b>
Select the input file (click "Choose * * * File")
Execute the operation (click "Import")

<b>1.1 Import Data ->Genotypic Data</b>
Choose one of different phenotypic formats (VCF, Fasta, P_link, Hapmap)

<b>1.2 Import Data ->Annotation</b>
Two annotation formats (GFF3 and BED) are supported,
The default format is GFF3. To select BED format annotation information, switch the radio box to BED

<b>1.3 Import Data ->Phenotypic Data</b>
Get phenotype file (Excel/WPS Save As ->Text file (tab separated) or CSV (comma separated))
Right radio box switch file format
When the lower check box is not checked (Show pheno data after import), only the header will be displayed after importing the file
Check the lower check box to output the file content to the operation interface

<b>1.4 Import Data ->Genotypic data</b>
At most two information files can be imported
It is recommended to integrate information other than phenotype into one file
Right radio box switch file format
When the lower check box is not checked (Show access information data after import), only the header will be displayed after importing the file
Check the lower check box to output the file content to the operation interface

<h3>2. Genotype pretreatment (optional) (Pre Process ->Filter Genotypic data)</h3>
There are four filtering modes:
  1. Filter by location
  2. Filter by type
  3. Filter by location and type
  4. No filtering
When filtering by position, you must specify three parameters: Chromosome, Start and End
When filtering by type, check at least one checkbox to specify the type of DNA segment to be retained
When not filtering, various parameters need not be adjusted

<h3>3. Haplotype identification and result preservation</h3>
<b>3.1 Haplotype ->Haplotype Identificition</b>
Select genotype file type (default is vcf)
Click "Identify Haplotye"
Click "Summary"
<b>3.2 Haplotype ->Save Result</b>
Click "Choose Result Directory" and select the file saving directory
Input file name
Select Save File Type
Click "Save" to save the result

<h3>4. Visualization</h3>
4.0 Basic steps:
 Adjust image size (unit: pixel) and resolution
 Adjust Optional Parameters
 Click (Plot IMG)
4.1 Overview of haplotype variation information (Visualization ->Haplotype Table)
4.2 Visualization ->Display Variables on Gene Model
4.3 haploid network (Visualization ->HaploNet)
4.4 Geographical distribution of main haplotypes (Visualization ->Geo distribution)
4.5 Phenotype Comparison
4.6 LD Heat Map (Visualization ->LD Heat Map)'
    h_ch <-"<h2>\u64cd\u4f5c\u6b65\u9aa4\uff1a</h2>
<h3>1. \u5bfc\u5165\u6570\u636e</h3>
<p>1.0 \u57fa\u672c\u6b65\u9aa4\uff1a
  \u9009\u62e9\u8f93\u5165\u6587\u4ef6\uff08\u70b9\u51fb\u201cChoose *** File\u201d\uff09
  \u6267\u884c\u64cd\u4f5c\uff08\u70b9\u51fb\u201cImport\u201d\uff09

<b>1.1 \u5bfc\u5165\u57fa\u56e0\u578b\u6587\u4ef6\uff08Import Data->Genotypic Data\uff09</b>
  \u4e0d\u540c\u8868\u578b\u683c\u5f0f\u4efb\u9009\u5176\u4e00\uff08VCF\uff0cFasta\uff0cP_link\uff0cHapmap\uff09

<b>1.2 \u5bfc\u5165\u6ce8\u91ca\u6587\u4ef6\uff08Import Data->Annotation\uff09</b>
  \u652f\u6301\u4e24\u79cd\u6ce8\u91ca\u683c\u5f0f\uff08GFF3\u548cBED\uff09\uff0c
  \u9ed8\u8ba4\u4e3aGFF3\u683c\u5f0f\uff0c\u5982\u9700\u9009\u62e9BED\u683c\u5f0f\u6ce8\u91ca\u4fe1\u606f\uff0c\u5c06\u5355\u9009\u6846\u5207\u6362\u4e3aBED\u5373\u53ef

<b>1.3 \u5bfc\u5165\u8868\u578b\u6587\u4ef6\uff08Import Data->Phenotypic Data\uff09</b>
  \u8868\u578b\u6587\u4ef6\u83b7\u5f97\uff08Excel/WPS \u53e6\u5b58\u4e3a -> \u6587\u672c\u6587\u4ef6\uff08\u5236\u8868\u7b26\u5206\u9694\uff09\u6216CSV\uff08\u9017\u53f7\u5206\u9694\uff09\uff09
  \u53f3\u4fa7\u5355\u9009\u6846\u5207\u6362\u6587\u4ef6\u683c\u5f0f
  \u4e0d\u52fe\u9009\u4e0b\u4fa7\u590d\u9009\u6846\u65f6\uff08Show pheno data after import\uff09\uff0c\u5bfc\u5165\u6587\u4ef6\u540e\u53ea\u5c55\u793a\u8868\u5934
  \u52fe\u9009\u4e0b\u4fa7\u590d\u9009\u6846\u53ef\u4ee5\u5c06\u6587\u4ef6\u5185\u5bb9\u8f93\u51fa\u5230\u64cd\u4f5c\u754c\u9762

<b>1.4 \u5bfc\u5165\u4e2a\u4f53\uff08Accession\uff09\u4fe1\u606f\uff08Import Data->Genotypic data\uff09</b>
  \u6700\u591a\u652f\u6301\u5bfc\u5165\u4e24\u4e2a\u4fe1\u606f\u6587\u4ef6
  \u63a8\u8350\u5c06\u9664\u8868\u578b\u5916\u7684\u5176\u4ed6\u4fe1\u606f\u6574\u5408\u5230\u4e00\u4e2a\u6587\u4ef6\u4e2d
  \u53f3\u4fa7\u5355\u9009\u6846\u5207\u6362\u6587\u4ef6\u683c\u5f0f
  \u4e0d\u52fe\u9009\u4e0b\u4fa7\u590d\u9009\u6846\u65f6\uff08Show accession information data after import\uff09\uff0c\u5bfc\u5165\u6587\u4ef6\u540e\u53ea\u5c55\u793a\u8868\u5934
  \u52fe\u9009\u4e0b\u4fa7\u590d\u9009\u6846\u53ef\u4ee5\u5c06\u6587\u4ef6\u5185\u5bb9\u8f93\u51fa\u5230\u64cd\u4f5c\u754c\u9762</p>
<h3>2. \u57fa\u56e0\u578b\u9884\u5904\u7406\uff08\u53ef\u9009\uff09\uff08Pre-Process -> Filter Genotypic data\uff09</h3>
<p>  \u5171\u56db\u79cd\u7b5b\u9009\u6a21\u5f0f\uff1a
    1) \u6309\u4f4d\u7f6e\u7b5b\u9009
    2) \u6309\u7c7b\u578b\u7b5b\u9009
    3) \u6309\u4f4d\u7f6e\u548c\u7c7b\u578b\u7b5b\u9009
    4) \u4e0d\u7b5b\u9009
  \u6309\u4f4d\u7f6e\u7b5b\u9009\u65f6\uff0c\u5fc5\u987b\u6307\u5b9a\u67d3\u8272\u4f53\uff08Chromosome\uff09\u3001\u5f00\u59cb\u4f4d\u7f6e\uff08Start\uff09\u3001\u7ec8\u6b62\u4f4d\u7f6e\uff08End\uff09\u4e09\u4e2a\u53c2\u6570
  \u6309\u7c7b\u578b\u7b5b\u9009\u65f6\uff0c\u81f3\u5c11\u52fe\u9009\u4e00\u4e2a\u590d\u9009\u6846\u6765\u6307\u5b9a\u4fdd\u7559\u7684DNA\u7247\u6bb5\u7c7b\u578b
  \u4e0d\u7b5b\u9009\u65f6\uff0c\u5404\u7c7b\u53c2\u6570\u4e0d\u7528\u8c03\u6574
</p>
<h3>3. \u5355\u500d\u578b\u9274\u5b9a\u53ca\u7ed3\u679c\u4fdd\u5b58</h3>
<p ><b>3.1 \u5355\u500d\u578b\u9274\u5b9a\uff08Haplotype->Haplotype Identificition\uff09</b>
  \u9009\u62e9\u57fa\u56e0\u578b\u6587\u4ef6\u7c7b\u578b\uff08\u9ed8\u8ba4\u4e3avcf\uff09
  \u70b9\u51fb\u201cIdentify Haplotye\u201d
  \u70b9\u51fb\u201cSummary\u201d

<b>3.2 \u5355\u500d\u578b\u7ed3\u679c\u4fdd\u5b58\uff08Haplotype->Save Result\uff09</b>
  \u70b9\u51fb\u201cChoose Result Directory\u201d\uff0c\u9009\u62e9\u6587\u4ef6\u4fdd\u5b58\u76ee\u5f55
  \u8f93\u5165\u6587\u4ef6\u540d\u79f0
  \u9009\u62e9\u4fdd\u5b58\u6587\u4ef6\u7c7b\u578b
  \u70b9\u51fb\u201cSave\u201d\uff0c\u4fdd\u5b58\u7ed3\u679c
</p>
<h3>4. \u53ef\u89c6\u5316\u4f5c\u56fe</h3>
<p><b>4.0 \u57fa\u672c\u6b65\u9aa4\uff1a</b>
 \u8c03\u6574\u56fe\u50cf\u5927\u5c0f\uff08\u5355\u4f4d\uff1a\u50cf\u7d20\uff09\u3001\u56fe\u50cf\u5206\u8fa8\u7387\uff08resolution\uff09
 \u8c03\u6574\u53ef\u9009\u53c2\u6570
 \u70b9\u51fb\uff08Plot IMG\uff09

<b>4.1 \u5355\u500d\u578b\u53d8\u5f02\u4fe1\u606f\u603b\u89c8\uff08Visualization->Haplotype Table\uff09</b>

<b>4.2 \u5355\u500d\u578b\u53d8\u5f02\u4f4d\u70b9\u5728\u57fa\u56e0\u4e0a\u7684\u4f4d\u7f6e\uff08Visualization->Display Variants on Gene Model\uff09</b>

<b>4.3 \u5355\u500d\u578b\u7f51\u7edc\uff08Visualization->HaploNet\uff09</b>

<b>4.4 \u4e3b\u8981\u5355\u500d\u578b\u5730\u7406\u5206\u5e03\uff08Visualization->Geo-distribution\uff09</b>

<b>4.5 \u8868\u578b\u6bd4\u8f83\uff08Visualization->Phenotype Comparison\uff09</b>

<b>4.6 LD\u70ed\u56fe\uff08Visualization->LD Heatmap\uff09</b>

</p>"
    h_ch <- stringr::str_replace_all(h_ch,"\n","<br/>")
    h_en <- stringr::str_replace_all(h_en,"\n","<br/>")

    ## server
    server <- function(input, output, session) {
        # Import Genotypic data
        # IMPORT_FASTA
        vcf <- seq <- p_link <- geno.gt  <- pheno <- gff <-
            accinfo <- hapResult <- hapSummary <- NULL

        # import vcf
        vcf <- eventReactive(input$vcf_im, {
            vcf <- if(is.surfix(input$vcf_in, c("vcf","vcf.gz")))
                geneHapR::import_vcf(input$vcf_in) else
                    NULL
            if(!is.null(vcf)){
                chrs <- as.list(unique(vcf@fix[,1]))
                names(chrs) <- as.vector(chrs)
                updateSelectInput(session, "filter_chr",
                                  choices = chrs,
                                  selected = chrs[[1]])
                start <- min(as.numeric(vcf@fix[,2]))
                end <- max(as.numeric(vcf@fix[,2]))
                updateNumericInput(session,"filter_start",
                                   value = start,
                                   min = start,
                                   max = end)
                updateNumericInput(session,"filter_end",
                                   value = end,
                                   min = start,
                                   max = end)
            }
            vcf
        }
        )
        observeEvent(input$vcf_bt,{
            file <- utils::choose.files(multi = FALSE, filters = ffilter[c("all", "vcf"),])
            updateTextInput(session, "vcf_in", value = file)
        })
        observeEvent(input$vcf_im,{
            info <- vcf()
            output$vcf_info <- renderPrint({info})
        })

        # import fasta
        fasta <- eventReactive(input$seq_im, {
            if(is.surfix(input$seq_in, c("fasta", "fa")))
                geneHapR::import_seqs(input$seq_in) else
                    NULL
        })
        # choose fatsa file
        observeEvent(input$seq_bt,{
            file <- utils::choose.files(multi = FALSE, filters = ffilter[c("all", "fa"),])
            updateTextInput(session, "seq_in", value = file)
        })
        # import fasta file
        observeEvent(input$seq_im,{
            info <- data.frame("seq_name" = names(fasta()),
                               "seq_length" = Biostrings::nchar(fasta()))
            output$seq_data <- renderTable({info})
        })


        # import p.link
        plink <- eventReactive(input$plink_im, {
            p.link <- if(is.surfix(input$map_in, "map") & is.surfix(input$ped_in, "ped") )
                geneHapR::import_plink.pedmap(pedfile = input$ped_in,
                                              mapfile = input$map_in,
                                              sep_ped = input$ped_sep,
                                              sep_map = input$map_sep) else
                                                  NULL
            if(!is.null(p.link)){
                chrs <- as.list(unique(p.link$map[,1]))
                names(chrs) <- as.vector(chrs)
                updateSelectInput(session, "filter_chr",
                                  choices = chrs,
                                  selected = chrs[[1]])
                start <- min(as.numeric(p.link$map[,4]))
                end <- max(as.numeric(p.link$map[,4]))
                updateNumericInput(session,"filter_start",
                                   value = start,
                                   min = start,
                                   max = end)
                updateNumericInput(session,"filter_end",
                                   value = end,
                                   min = start,
                                   max = end)
            }
            p.link
        })
        # choose ped file
        observeEvent(input$ped_bt, {
            file <- utils::choose.files(filters = ffilter[c("all", "plink"),], multi = FALSE)
            updateTextInput(session, "ped_in", value = file)
        })
        # choose map file
        observeEvent(input$map_bt, {
            file <- utils::choose.files(filters = ffilter[c("all", "plink"),], multi = FALSE)
            updateTextInput(session, "map_in", value = file)
        })
        # import plink files
        observeEvent(input$plink_im, {
            plink_object <- plink()
            output$plink_info <- renderPrint({
                cat("Indeividuals: ", nrow(plink_object$ped),"\n")
                cat("Markers: ", nrow(plink_object$map))
            })
        })

        # import annotation files
        gff <- eventReactive(input$gff_im,{
            gff <- if(is.surfix(input$gff_in, c("gff","gff3","bed","bed6","bed4")))
                switch(input$gff_format,
                       "GFF" = geneHapR::import_gff(input$gff_in),
                       "BED" = geneHapR::import_bed(input$gff_in)) else NULL
            if(!is.null(gff)){
                choices <- as.vector(unique(gff$type))
                updateCheckboxGroupInput(session,
                                         "filter_type",
                                         choiceNames = choices,
                                         choiceValues = choices)
                updateCheckboxGroupInput(session,
                                         "geneElement", inline = TRUE,
                                         choiceNames = choices,
                                         choiceValues = choices)

            }
            gff
        })
        # choose annotation file
        observeEvent(input$gff_bt,{
            file <- utils::choose.files(multi = FALSE, filters = ffilter[c("all", "bed","gff"),])
            updateTextInput(session, "gff_in", value = file)
        })
        # import annotation file
        observeEvent(input$gff_im,{
            output$annotation_data <- renderPrint(gff())
        })

        pheno <- eventReactive(input$pheno_im,{
            if(any(is.surfix(input$pheno_in, c("txt","csv")))){
                switch(input$pheno_format,
                       "txt" = geneHapR::import_AccINFO(input$pheno_in),
                       "csv" = utils::read.csv(input$pheno_in, header = TRUE, row.names = 1))
            }
        })
        observeEvent(input$pheno_bt,{
            file <- utils::choose.files(multi = FALSE, filters = ffilter[c("all","csv","txt"),])
            updateTextInput(session, "pheno_in", value = file)
        })
        observeEvent(input$pheno_im,{
            pheno_object <- pheno()
            output$pheno_info <- renderUI({
                text <- paste(h4("Imported pheno names:"), br(),
                              paste(names(pheno_object), collapse = "; "))
                HTML(text)})
            if(input$pheno_showdata)
                output$pheno_data <- renderTable(pheno_object) else
                    output$pheno_data <- renderTable(NULL)
            choices <- as.list(names(pheno_object))
            names(choices) <- names(pheno_object)
            updateSelectInput(session, "hapvspheno_phenoname", choices = choices)
        })


        accinfo <- eventReactive(input$accinfo_im,{
            if(nchar(input$accinfo_in) > 0 & any(is.surfix(input$accinfo_in, c("txt","csv")))){
                switch(input$accinfo_format,
                       "txt" = geneHapR::import_AccINFO(input$accinfo_in),
                       "csv" = utils::read.csv(input$accinfo_in, header = TRUE, row.names = 1))
            } else NULL
        })
        observeEvent(input$accinfo_bt,{
            file <- utils::choose.files(multi = FALSE, filters = ffilter[c("all","csv","txt"),])
            updateTextInput(session, "accinfo_in", value = file)
        })
        observeEvent(input$accinfo_im,{
            accinfo_object <- accinfo()
            output$accinfo_info <- renderUI({
                text <- paste(h4("Names: "), br(),
                              paste(names(accinfo_object), collapse = "; "))
                HTML(text)})
            if(input$accinfo_showdata)
                output$accinfo_data <- renderTable(head(accinfo_object)) else
                    output$accinfo_data <- renderTable(NULL)
            choices <- as.list(names(accinfo()), names(accinfo1()))
            choices <- c(choices, "none" = "none")
            names(choices) <- as.vector(choices)
            updateSelectInput(session, "plothapnet_group",
                              choices = choices, selected = "none")
            updateSelectInput(session, "hapdistribution_loncol",
                              choices = choices, selected = "none")
            updateSelectInput(session, "hapdistribution_latcol",
                              choices = choices, selected = "none")
        })
        accinfo1 <- eventReactive(input$accinfo_im1,{
            if(any(is.surfix(input$accinfo_in1, c("txt","csv")))){
                switch(input$accinfo_format1,
                       "txt" = geneHapR::import_AccINFO(input$accinfo_in1),
                       "csv" = utils::read.csv(input$accinfo_in1, header = TRUE, row.names = 1))
            } else NULL
        })
        observeEvent(input$accinfo_bt1,{
            file <- utils::choose.files(multi = FALSE, filters = ffilter[c("all","csv","txt"),])
            updateTextInput(session, "accinfo_in1", value = file)
        })
        observeEvent(input$accinfo_im1,{
            accinfo_object <- accinfo1()
            output$accinfo_info1 <- renderUI({
                text <- paste(h4("Names: "), br(),
                              paste(names(accinfo_object), collapse = "; "))
                HTML(text)})
            if(input$accinfo_showdata)
                output$accinfo_data1 <- renderTable(head(accinfo_object)) else
                    output$accinfo_data1 <- renderTable(NULL)
            choices <- as.list(names(accinfo()), names(accinfo1()))
            choices <- c(choices, "none" = "none")
            names(choices) <- as.vector(choices)
            updateSelectInput(session, "plothapnet_group",
                              choices = choices, selected = "none")
            updateSelectInput(session, "hapdistribution_loncol",
                              choices = choices, selected = "none")
            updateSelectInput(session, "hapdistribution_latcol",
                              choices = choices, selected = "none")
        })


        fvcf <- eventReactive(input$hapresult_bt,{
            vcf <- vcf()
            message("orignal variants: ", nrow(vcf@fix))
            message("filter mode: ", input$filter_mode)
            if(input$filter_mode == "none") vcf <- vcf() else {
                if(input$filter_mode %in% c("POS","both")) {
                    message("filter position: ", input$filter_chr,": ",
                            input$filter_start,"-",
                            input$filter_end)
                }
                if(input$filter_mode %in% c("type","both")) {
                    message("filter type: ", input$filter_type)
                    if(is.null(gff()) | is.null(input$filter_type))
                        output$hapresult <- renderPrint("Please input annotation or cancel filter by type")
                }
                vcf <- geneHapR::filter_vcf(vcf(), gff = gff(),
                                            mode = input$filter_mode,
                                            Chr = input$filter_chr,
                                            start = input$filter_start,
                                            end = input$filter_end,
                                            type = input$filter_type)
                message("variants after filter: ", nrow(vcf@fix))
            }
            vcf
        })

        fplink <- eventReactive(input$hapresult_bt,{
            plink <- plink()
            message("orignal variants: ", nrow(plink$map))
            message("filter mode: ", input$filter_mode)

            if(input$filter_mode == "none") plink() else{
                if(input$filter_mode %in% c("POS","both")) {
                    message("filter position: ", input$filter_chr,": ",
                            input$filter_start,"-",
                            input$filter_end)
                }
                if(input$filter_mode %in% c("type","both")) {
                    message("filter type: ", input$filter_type)
                    if(is.null(gff()) | is.null(input$filter_type))
                        output$hapresult <- renderPrint("Please input annotation or cancel filter by type")
                }

                plink <- geneHapR::filter_plink.pedmap(plink(), gff = gff(),
                                                       mode = input$filter_mode,
                                                       Chr = input$filter_chr,
                                                       start = input$filter_start,
                                                       end = input$filter_end,
                                                       type = input$filter_type)

                message("variants after filter: ", nrow(plink$map))
            }
            plink
        })


        # haplotyping
        geno_fmt <- reactive({input$hapresult_source})
        hapre <- eventReactive(input$hapresult_bt,{
            hapre <- try(switch(geno_fmt(),
                                "vcf" = geneHapR::vcf2hap(fvcf(),
                                                          hapPrefix = input$happrefix,
                                                          hyb_remove = input$hyb_remove,
                                                          na.drop = input$na_drop),
                                "Fasta" = geneHapR::seqs2hap(fasta(),
                                                             hapPrefix = input$happrefix,
                                                             hyb_remove = input$hyb_remove,
                                                             na.drop = input$na_drop),
                                "p.link" = geneHapR::plink.pedmap2hap(fplink(),
                                                                      hapPrefix = input$happrefix,
                                                                      hyb_remove = input$hyb_remove,
                                                                      na.drop = input$na_drop),
                                "Table" = NULL,
                                "hapmap" = NULL))
            hapre
        })
        hapsum <- eventReactive(input$hapsummary_bt,{
            geneHapR::hap_summary(hapre())

        })


        observeEvent(input$hapresult_bt,{
            hapresult_object <- hapre()
            output$hapresult <- renderPrint(hapresult_object)
            # output$hapsummary <- renderPrint(hapsummary_object)
            if(inherits(hapresult_object, "hapResult")){
                haps <- names(attr(hapresult_object, "freq"))
                updateRadioButtons(session, inputId = "hapvspheno_hap1",
                                   choiceNames = haps,choiceValues = haps)
                updateRadioButtons(session, inputId = "hapvspheno_hap2",
                                   choiceNames = haps,choiceValues = haps)
                haps <- haps[attr(hapresult_object, "freq") >= 2]
                updateCheckboxGroupInput(session, inputId = "hapdistribution_hapnames",
                                         choiceNames = haps,choiceValues = haps, inline = TRUE)
                info <- hapresult_object[3, -c(1, ncol(hapresult_object))] %>%
                    data.frame() %>% t()
                info <- info[,1]
                info <- sapply(info, function(x) strsplit(x, "=")[[1]][1]) %>%
                    unlist() %>% unique()
                if(length(info) >= 1){
                    infos <- as.list(c(info,"none"))
                    names(infos) <- c(info, "none")
                    updateSelectInput(session, "plothaptable_infotag",
                                      choices = infos,
                                      selected = "none")
                }
            }
        })

        observeEvent(input$hapsummary_bt,{
            hapsummary_object <- hapsum()
            output$hapsummary <- renderPrint(hapsummary_object)

            chr <- unique(hapsummary_object[1, c(2:(ncol(hapsummary_object) - 2))])
            chr_list <- as.list(chr)
            names(chr_list) <- chr
            updateSelectInput(session, "displayVarOnGeneModel_chromosome",
                              choices = chr_list)
            updateSelectInput(session, "LD_heatmap_Chr",
                              choices = chr_list)
            pos <- as.numeric(hapsummary_object[2,c(2:(ncol(hapsummary_object) - 2))])
            updateNumericInput(session, "displayVarOnGeneModel_start", value = min(pos))
            updateNumericInput(session, "displayVarOnGeneModel_end", value = max(pos))
            updateNumericInput(session, "LD_heatmap_start", value = min(pos))
            updateNumericInput(session, "LD_heatmap_end", value = max(pos))
        })

        # save haplotype result
        observeEvent(input$chooseresultdir,{
            dir <- choose.dir(default = getwd())
            updateTextInput(session, "resultdir", value = dir)
        })
        observeEvent(input$download_hapresult, {
            resdir <- input$resultdir
            resultname <- tolower(input$resultname)
            filepath <- switch (input$saveresult,
                                "hapSummary" = {
                                    if(is.surfix(resultname,c("csv", "txt","hapsummary")))
                                        paste0(resdir, "/", input$resultname) else
                                            paste0(resdir, "/", input$resultname, ".hapSummary")
                                },
                                "hapResult" = {
                                    if(is.surfix(resultname,c("csv", "txt","hapresult")))
                                        paste0(resdir, "/", input$resultname) else
                                            paste0(resdir, "/", input$resultname, ".hapResult")
                                })
            object <- switch(input$saveresult,
                             "hapSummary" = hapsum(),
                             "hapResult" = hapre())

            geneHapR::write.hap(object, file = filepath)

        })

        # Visulization
        # plotHapTable
        observeEvent(input$plothaptable,{
            genename <- if(input$plothaptable_genename == "") FALSE else
                input$plothaptable_genename
            title <- input$plothaptable_title
            if(nchar(title) == 0) title <- ""
            if(input$plothaptable_infotag == "none") {
                infotag <- tagsplit <- tagfield <- tagname <- NULL
            } else {
                infotag <- input$plothaptable_infotag
                tagsplit <- input$plothaptable_tagsplit
                tagfield <- input$plothaptable_tagfield
                tagname <- input$plothaptable_tagname
            }
            hapSummary <- hapsum()
            output$plothaptable <- renderPlot(
                {
                    geneHapR::plotHapTable(
                        hapSummary = hapSummary,
                        title = title,
                        hapPrefix = input$plothaptable_prefix,
                        geneName = genename,
                        INFO_tag = infotag,
                        tag_split = tagsplit,
                        tag_field = tagfield,
                        tag_name = tagname,
                        angle = input$plothaptable_angle)
                },
                width = input$img_width,
                height = input$img_height,
                res = input$img_res)})

        # plothapnet
        legend_hapnet <- reactive(
            {if(any(input$plothapnet_legend_size,
                    input$plothapnet_legend_color))
                c(input$plothapnet_x, input$plothapnet_y) else
                    FALSE
            })
        observeEvent(input$plothaplonet,{
            message("start")

            # accinfo1 <- accinfo1()
            # accinfo <- accinfo1()
            # message("startprocessing")
            # if(input$plothapnet_group %in% names(accinfo1))
            #   accinfo <- accinfo1
            # message("processing")
            if(input$plothapnet_group != "none"){
                accinfo <- accinfo()
                message("accinfo")
                hapNet <- geneHapR::get_hapNet(hapsum(), AccINFO = accinfo,
                                               groupName = input$plothapnet_group)
            } else {
                hapNet <- geneHapR::get_hapNet(hapsum())
            }
            message("ploting")

            # message("lenged: ", legend)
            xlim <- input$plothapnet_xlim
            ylim <- input$plothapnet_ylim
            f <- tempfile()
            png(f)
            geneHapR::plotHapNet(hapNet)
            m <- round(par("usr"))
            dev.off()
            unlink(f)
            # updateSliderInput(session, "plothapnet_xlim",
            #                   min = m[1]*3,
            #                   max = m[2]*3,
            #                   value = c(m[1],m[2]))
            # updateSliderInput(session, "plothapnet_ylim",
            #                   min = m[3]*3,
            #                   max = m[4]*3,
            #                   value = c(m[3],m[4]))
            # updateSliderInput(session, "plothapnet_x",
            #                   min = m[1]*3,
            #                   max = m[2]*3,
            #                   value = m[1] + m[2])
            # updateSliderInput(session, "plothapnet_y",
            #                   min = m[3]*3,
            #                   max = m[4]*3,
            #                   value = m[3] + m[4])



            output$plothaplonet <- renderPlot({
                geneHapR::plotHapNet(hapNet, scale = input$plothapnet_scale,
                                     show.mutation = input$plothapnet_showmutation,
                                     cex = input$plothapnet_cex,
                                     label = input$plothapnet_labels,
                                     cex.legend = input$plothapnet_cexlegend,
                                     col.link = input$plothapnet_collink,
                                     link.width = input$plothapnat_linkwidth,
                                     legend = legend_hapnet(),
                                     show_color_legend = input$plothapnet_legend_color,
                                     show_size_legend = input$plothapnet_legend_size,
                                     xlim = input$plothapnet_xlim,
                                     ylim = input$plothapnet_ylim)
            }, width = input$img_width, height = input$img_height, res = input$img_res)


            message("end")

        })

        # display vars on model
        observeEvent(input$displayvar, {
            message("start ploting vars")
            Chr <- input$displayVarOnGeneModel_chromosome
            startPOS <- input$displayVarOnGeneModel_start
            endPOS <- input$displayVarOnGeneModel_end
            type <- input$displayVarOnGeneModel_type
            cex <- input$displayVarOnGeneModel_cex
            geneElement <- input$geneElement
            gff <- gff()
            hapSummary <-  hapsum()
            output$displayvarongenemodel <- renderPlot({
                geneHapR::displayVarOnGeneModel(gff = gff,hapSummary =  hapSummary,
                                                Chr = Chr,
                                                startPOS = startPOS,
                                                endPOS = endPOS,
                                                type = type,
                                                cex = cex,
                                                geneElement = geneElement)
            }, width = input$img_width, height = input$img_height, res = input$img_res)
            message("ending plot vars")
        })

        # hapdistribution
        observeEvent(input$hapdistribution_legend,{
            if(input$hapdistribution_legend == "other"){
                show("hapdistribution_x")
                show("hapdistribution_y")
            } else {
                hide("hapdistribution_x")
                hide("hapdistribution_y")
            }
        })
        observeEvent(input$hapdistribution_region,{
            if(input$hapdistribution_region == "other"){
                show("hapdistribution_region_other")
            } else {
                hide("hapdistribution_region_other")
            }
        })
        observeEvent(input$hapdistribution,{
            message("begin")
            # library(mapdata)
            # library(maptools)
            # accinfo1 <- accinfo1()
            # if(input$hapdistribution_loncol %in% names(accinfo1))
            #   accinfo <- accinfo1 else
            accinfo <- accinfo()
            if(input$hapdistribution_legend == "other")
                legend <- c(input$hapdistribution_x, input$hapdistribution_y) else
                    legend <- input$hapdistribution_legend
            if(input$hapdistribution_region == "other")
                database <- input$hapdistribution_region_other else
                    database <- input$hapdistribution_region
            message("ploting")
            LON.col = input$hapdistribution_loncol
            LAT.col = input$hapdistribution_latcol

            hap <- hapre()
            output$hapdistribution <- renderPlot({
                if(input$hapdistribution_latcol == "none" ||
                   input$hapdistribution_loncol == "none") {
                    plot.new()
                    title("Please choose latitude and longitude")
                } else
                    geneHapR::hapDistribution(hap = hap, AccINFO = accinfo,
                                              LON.col = LON.col,
                                              LAT.col = LAT.col,
                                              hapNames = input$hapdistribution_hapnames,
                                              database = database,
                                              regions = ".",
                                              legend = legend,
                                              lwd.pie = input$hapdistribution_lwdpie,
                                              lwd = input$hapdistribution_lwd,
                                              cex.legend = input$hapdistribution_cexlegend,
                                              symbolSize = input$hapdistribution_symbolsize,
                                              xlim = input$hapdistribution_xlim,
                                              ylim = input$hapdistribution_ylim)
            }, width = input$img_width, height = input$img_height, res = input$img_res)
            message("end")

        })

        # hapVsPhenos
        observeEvent(input$hapvaspheno,{
            result <- try(geneHapR::hapVsPheno(hap = hapre(), pheno = pheno(),
                                               mergeFigs = FALSE,
                                               minAcc = input$hapvspheno_minacc,
                                               phenoName = input$hapvspheno_phenoname,
                                               title = input$hapvspheno_title,
                                               hapPrefix = input$hapvspheno_prefix,
                                               angle = input$hapvspheno_angle,
                                               outlier.rm = input$hapvspheno_outlierm))
            if(inherits(result, "try-error")){
                output$hapvaspheno <- renderPlot({
                    plot.new()
                    title(as.character(result))
                })
            } else
                output$hapvaspheno <- renderPlot({result$fig_Violin})

        })

        # LD_heatmap
        observeEvent(input$LD_heatmap_addmap,{
            if(input$LD_heatmap_addmap){
                show(id = "LD_heatmap_Chr")
                show(id = "LD_heatmap_start")
                show(id = "LD_heatmap_end")
                show(id = "LD_heatmap_geneID")
                show(id = "LD_heatmap_maplocation")
                show(id = "LD_heatmap_mapheight")
            } else {
                hide(id = "LD_heatmap_Chr")
                hide(id = "LD_heatmap_start")
                hide(id = "LD_heatmap_end")
                hide(id = "LD_heatmap_geneID")
                hide(id = "LD_heatmap_maplocation")
                hide(id = "LD_heatmap_mapheight")
            }
        })
        observeEvent(input$plotLD,{
            title <- if(nchar(input$LD_heatmap_title) == 0) " " else
                input$LD_heatmap_title
            LDmeasure <- input$LD_heatmap_LDmeasure
            distances <- input$LD_heatmap_distance
            colorLegend <- input$LD_heatmap_colorLegend
            addmap <- input$LD_heatmap_addmap
            text <- input$LD_heatmap_text

            if(addmap){
                Chr <- input$LD_heatmap_Chr
                start <- input$LD_heatmap_start
                end <- input$LD_heatmap_end
                color <- input$LD_heatmap_clolor
                if(color == "grey") color <- grDevices::grey.colors(40) else color <- NULL
                geneMapLocation <- input$LD_heatmap_maplocation
                map.height <- input$LD_heatmap_mapheight
                hapre <- hapre()
                gff <- gff()
                output$plotLD <- renderPlot({
                    plot_LDheatmap(hapre, gff = gff, Chr = Chr,
                                   start = start, end = end,
                                   color = color,
                                   distances = distances,
                                   LDmeasure = LDmeasure,
                                   title = title,
                                   add.map = addmap,
                                   map.height = map.height,
                                   geneMapLocation = geneMapLocation,
                                   colorLegend = colorLegend,
                                   text = text)
                })
            } else {
                color <- input$LD_heatmap_clolor
                if(color == "grey") color <- grDevices::grey.colors(40) else color <- NULL
                hapre <- hapre()
                output$plotLD <- renderPlot({
                    plot_LDheatmap(hapre,
                                   color = color,
                                   distances = distances,
                                   LDmeasure = LDmeasure,
                                   title = title,
                                   add.map = addmap,
                                   colorLegend = colorLegend,
                                   text = text)
                })
            }
        })

        # working directory
        observeEvent(input$changewd,{
            setwd(choose.dir(default = getwd()))
        })
        # Language
        observeEvent(input$changelanguage,{
            if(input$changelanguage %% 2 == 1){
                updateActionButton(session, "changelanguage",
                                   label = "\u5207\u6362\u8BED\u8A00")
            } else {
                updateActionButton(session, "changelanguage",
                                   label = "Change Language")
            }
            toggle("help_en")
            toggle("help_ch")
        })


        root <- c(wd=getwd(),"root"="/")
        # shinyFiles::shinyDirChoose(input, "shinychooseresultdir", root = root,FALSE)
        observeEvent(input$shinychooseresultdir, {
            if(is.list(input$shinychooseresultdir)){
                updateTextInput(session, "resultdir", value = input$shinychooseresultdir)
                p <- input$shinychooseresultdir
                print(input$shinychooseresultdir)
                message(class(input$shinychooseresultdir))
                print(names(p))
                print(p$path)
                message(class(p))
            }
        })
    }



    ui <- function(...){

        tabPanel_welcome <- tabPanel( # Welcome Page
            "Welcompage",
            tabPanelBody(
                "Welcomepagenody",

                fluidRow(
                    column(width = 3,
                           actionButton("changewd", "Change Working Directory"),
                           actionButton("changelanguage", "Change Language")),
                    column(width = 7,
                           br(),
                           h1(id = "mainID", "Welcome to geneHapR-GUI"),
                           column(
                               width = 12, id = "help_en",
                               HTML(text = paste0("<a href='",
                                                  "https://gitee.com/zhangrenl/genehapr/wikis/An_example_in_Setaria_italica",
                                                  "'>","An Example in Setaria italica",
                                                  "</a>")),
                               p(HTML(text = h_en))),
                           hidden(column(
                               width = 12, id = "help_ch",
                               HTML(text = paste0("<a href='",
                                                  "https://gitee.com/zhangrenl/genehapr/wikis/An_example_in_Setaria_italica",
                                                  "'>","\u5355\u500d\u578b\u5206\u6790\u6848\u4f8b\uff08\u8c37\u5b50\uff09",
                                                  "</a>")),

                               p(HTML(text = h_ch))))
                           # HTML(stringr::str_replace_all(h_ch,"\n","<br/>"))
                    )
                )
            )
        )

        tabPanel_preprocess <- tabPanel(
            "Pre-Process",
            tabsetPanel(
                tabPanel(
                    "Filter Genotypic Data",
                    br(),
                    tabPanelBody(
                        "filter_vcf",
                        fluidRow(column(width = 3,
                                        radioButtons("filter_mode", label = "Filter variants by",
                                                     selected = "none",
                                                     choiceNames = c("Position", "Type",
                                                                     "Both of above","Do not filter"),
                                                     choiceValues = c("POS","type","both","none")),
                                        checkboxGroupInput("filter_type",
                                                           label = "Type",
                                                           choiceNames = c("CDS", "UTR"),
                                                           choiceValues = c("CDS", "UTR"),selected = "CDS")),
                                 column(width = 3,
                                        selectInput("filter_chr", label = "Chromosome",
                                                    choices = list()),
                                        numericInput("filter_start", label = "Start", value = 0),
                                        numericInput("filter_end", label = "End",value = 0))))),
                tabPanel("Extract Large VCF",
                         tabPanelBody(
                             "filter_large_vcf",
                             br(),
                             fluidRow(
                                 column(width = 2,
                                        actionButton("lvcf_filter_vcfin_bt", width = 130, label = "Input vcf file")),
                                 column(width = 6,
                                        textInput("lvcf_filter_vcfin", NULL))),
                             fluidRow(
                                 column(width = 2,
                                        actionButton("lvcf_filter_vcfout_bt",width = 130, label = "Output vcf file")),
                                 column(width = 6,
                                        textInput("lvcf_filter_vcfout", NULL))),
                             br(),
                             fluidRow(column(width = 2,
                                             textInput("lvcf_filter_chr", label = "Chromosome")),
                                      column(width = 4,
                                             textInput("lvcf_filter_end", label = "Start"),
                                             textInput("lvcf_filter_start", label = "End"))),
                             actionButton("lvcf_filter_bt", "Extract From Large VCF", align = "center"))))
        )


        Geo_distribution <- tabPanel(
            "Geo-distribution",
            br(),
            fluidRow(
                column(width = 6,
                       selectInput("hapdistribution_loncol","Column name of longitude",
                                   choices = list("none" = "none")),
                       numericInput("hapdistribution_symbolsize", "Synbol size",
                                    step = 0.1, value = 1, min = 0, max = 10),
                       numericInput("hapdistribution_lwd", "Line width of map", value = 1),
                       selectInput("hapdistribution_region", "Region", selected = "china",
                                   choices = list("china" = "china",
                                                  "world" = "world",
                                                  "other" = "other")),
                       hidden(
                           textInput("hapdistribution_region_other", label = NULL)),
                ),
                column(width = 6,
                       selectInput("hapdistribution_latcol","Column name of latitude",
                                   choices = list("none" = "none")),
                       numericInput("hapdistribution_lwdpie", "Line width of pie", value = 0.5),
                       numericInput("hapdistribution_cexlegend", "Size of legend", value = 0.8),
                       selectInput("hapdistribution_legend", "Legend position",
                                   choices = list("none" = FALSE, 'left' = 'left',
                                                  'right' = 'right', 'top' = 'top', 'bottom' ='bottom',
                                                  'topleft' = 'topleft', 'topright' = 'topright',
                                                  'bottomright' = 'bottmoright', 'bottomleft' = 'bottomleft',
                                                  "other" = "other")),
                       hidden(
                           numericInput("hapdistribution_x", "x", step = 1, width = 300, value = -170),
                           numericInput("hapdistribution_y", "y", step = 1, width = 200, value = -60)))),
            fluidRow(
                column(width = 6,
                       sliderInput("hapdistribution_xlim", "Longitude range",
                                   min = -180, max = 180, value =c(-180,180)),

                ),
                column(width = 6,
                       sliderInput("hapdistribution_ylim", "Latitude range",
                                   min = -90, max = 90, value =c(-90,90)))),
            checkboxGroupInput("hapdistribution_hapnames",
                               "Display haplotypes",inline = TRUE,
                               choiceNames = c("H001", "H002", "H003"),
                               choiceValues = c("H001", "H002", "H003"),
                               selected = c("H001", "H002", "H003")),
            actionButton("hapdistribution", "Plot IMG", width = "85%"),
            plotOutput("hapdistribution")
        )

        Phenotype_Comparison <- tabPanel(
            "Phenotype Comparison",
            br(),
            fluidRow(
                column(width = 6,
                       selectInput("hapvspheno_phenoname", "Pheno name",
                                   choices = list("first pheno" = "1", "second pheno" = 2),
                                   selected = "2"),
                       textInput("hapvspheno_prefix", "Prefix of haplotype names", value = "H"),
                       textInput("hapvspheno_title", "Title", value = "")),
                column(width = 6,
                       numericInput("hapvspheno_minacc", "Minimum number of accession", value = 5),
                       numericInput("hapvspheno_angle", "Angle of x labels", value = 0))),
            actionButton("hapvaspheno","Plot IMG", width = "85%"),
            checkboxInput("hapvspheno_outlierm", "Remove outlier", value = TRUE),
            #
            # h4("Custom comparisons"),
            # fluidRow(
            #   column(width = 6,
            #          actionButton("hapvspheno_addcompare", "Add"),
            #          radioButtons("hapvspheno_hap1", "",
            #                       choiceNames = c("H001","H002"),
            #                       choiceValues = c("H001","H002"))),
            #   column(width = 6,
            #          actionButton("hapvspheno_addcompare", "Clear"),
            #          radioButtons("hapvspheno_hap2", "",
            #                       choiceNames = c("H001","H002"),
            #                       choiceValues = c("H001","H002")))),
            br(),
            plotOutput("hapvaspheno"),
        )

        haploNet <- tabPanel(
            "Haplotype Network",
            br(),
            fluidRow(
                column(width = 4,
                       selectInput("plothapnet_scale", label = "Scale",selected = "log2",
                                   choices = list("none" = 1, "log10" = "log10", "log2" = "log2"))),
                column(width = 4,
                       selectInput("plothapnet_showmutation", label = "Mutation symbol",selected = "line",
                                   choices = list("none" = 0, "line" = "1",
                                                  "dot" = "2", "number of mutantions" = 3))),
                column(width = 4,
                       selectInput("plothapnet_group", "Group by", choices = list("none"="none")))),
            fluidRow(column(width = 6,
                            numericInput("plothapnet_collink", "Link color", value = 1)),
                     column(width = 6,
                            numericInput("plothapnat_linkwidth", "Line width", value = 1))),

            fluidRow(
                column(width = 6,
                       sliderInput("plothapnet_cex",
                                   label = "Size of label", min = 0, max = 3, value = 0.8, step = 0.05)),
                column(width = 6,
                       sliderInput("plothapnet_cexlegend",
                                   label = "Size of legend", min = 0, max = 3, value = 0.8, step = 0.05))),
            fluidRow(
                column(width = 6,
                       sliderInput("plothapnet_xlim", "X limit range",
                                   min = -200, max = 200, value =c(-1,1))
                ),
                column(width = 6,
                       sliderInput("plothapnet_ylim", "Y limit range",
                                   min = -200, max = 200, value =c(-1,1))
                )),
            h4("Legend options"),
            fluidRow(
                column(width = 6,
                       checkboxInput("plothapnet_legend_size", "Show size legend", value = TRUE)),
                column(width = 6,
                       checkboxInput("plothapnet_legend_color", "Show color legend", value = TRUE))),
            fluidRow(column(width = 6,
                            sliderInput("plothapnet_x", "x", min = -200, max = 200, value = -10)),
                     column(width = 6,
                            sliderInput("plothapnet_y", "y", min = -200, max = 200, value = -10))),
            br(),

            actionButton("plothaplonet", "Plot IMG", width = "85%"),
            checkboxInput("plothapnet_labels", "Display haplotype names", value = TRUE),

            br(),
            plotOutput("plothaplonet"),

        )

        Dispaly_Variants_on_Gene_Model <- tabPanel(
            "Dispaly Variants on Gene Model",
            br(),
            fluidRow(
                column(width = 6,
                       selectInput("displayVarOnGeneModel_chromosome","Chromosome",choices = list()),
                       numericInput("displayVarOnGeneModel_start", "Start", value = 0),
                       numericInput("displayVarOnGeneModel_end", "End", value = 0)),
                column(width = 6,
                       selectInput("displayVarOnGeneModel_type", label = "Variants type",
                                   choices = list("circle" = "circle", "pie" = "pie", "pin" = "pin", "flag" = "flag"),
                                   selected = "pin"),
                       sliderInput("displayVarOnGeneModel_cex",
                                   label = "Size of label", min = 0, max = 3, value = 0.8, step = 0.05),
                       checkboxGroupInput("geneElement",label = "Display elements", inline = TRUE,
                                          choiceNames = c("CDS","UTR"), choiceValues = c("CDS","UTR")))),
            br(),
            actionButton("displayvar", "Plot IMG", width = "85%"),
            plotOutput("displayvarongenemodel", height = "800px"),

        )

        Hapotype_Table <- tabPanel(
            "Hapotype Table",
            br(),
            fluidRow(
                column(width = 6,
                       textInput("plothaptable_prefix","Prefix of Haplotype Names",value = "H"),
                       textInput("plothaptable_title", "Title", value = " "),
                       textInput("plothaptable_genename", "Gene name", value = ""),
                       sliderInput("plothaptable_angle", label = "Angle of x-axis label",
                                   min = 0, max = 90, step = 45, value = 0)),


                column(width = 6,
                       selectInput("plothaptable_infotag", "Tag name in INFO", choices = list()),
                       textInput("plothaptable_tagname", "Tag name in IMG", value = ""),
                       textInput("plothaptable_tagsplit", "Tag split", value = "|"),
                       numericInput("plothaptable_tagfield", "Tag field", value = "1"))),
            actionButton("plothaptable", "Plot IMG", width = "85%"),
            br(),
            plotOutput("plothaptable", height = 750))

        plot_LDheatmap <- tabPanel(
            "LD Heatmap",
            br(),
            fluidRow(
                column(width = 6,
                       textInput("LD_heatmap_geneID", "Gene ID", " "),
                       selectInput("LD_heatmap_distance","Distance type",
                                   choices = list("physical" = "physical",
                                                  "genetic" = "genetic"),
                                   selected = "physical"),
                       textInput("LD_heatmap_title", "Title", value = " "),
                       sliderInput("LD_heatmap_mapheight", "Map height",
                                   min = 0.0025, max = 0.1,
                                   step = 0.0025,
                                   value = 0.02)),
                column(width = 6,
                       fluidRow(
                           column(
                               width = 4,
                               selectInput("LD_heatmap_Chr","Chromosome",
                                           choices = list("none" = "none"))
                           ),
                           column(
                               width = 4,
                               numericInput("LD_heatmap_start",
                                            "Start", step = 500,
                                            value = 1)),
                           column(
                               width = 4,
                               numericInput("LD_heatmap_end",
                                            "End", step = 500,
                                            value = 1)
                           )),
                       selectInput("LD_heatmap_clolor", "Color",
                                   choices = list("whiteyellowred"="whiteyellowred",
                                                  "grey" = "grey")),
                       selectInput("LD_heatmap_LDmeasure","LDmeasure",
                                   choices = list("allelic correlation r^2" = "r",
                                                  "Lewontin's |D'|" = "D'"),
                                   selected = "physical"),
                       sliderInput("LD_heatmap_maplocation", "Map location",
                                   min = 0, max = 1, value = 0.15, step = 0.025))),
            actionButton("plotLD", "Plot IMG", width = "85%"),
            checkboxInput("LD_heatmap_colorLegend", "Show color legend", TRUE),
            checkboxInput("LD_heatmap_addmap", "Show gene map", FALSE),
            checkboxInput("LD_heatmap_text", "Add text to cell", FALSE),

            br(),
            plotOutput("plotLD", height = 750))
        tabPanel_Visualization <- tabPanel( # visualization start
            "Visualization",
            # siderbar UI start
            sidebarPanel(
                width = 3,
                h3("Image  Options"),
                # textInput("img_name", "File name")
                numericInput("img_width", "width", value = 960,min = 1,step = 1),
                numericInput("img_height", "height", value = 720,min = 1,step = 1),
                numericInput("img_res", "resolution ", value = 72,min = 1,step = 1)),
            # siderbar UI End

            # plotHap table UI start
            mainPanel(tabsetPanel(
                Hapotype_Table,
                Dispaly_Variants_on_Gene_Model,
                haploNet,
                Geo_distribution,
                Phenotype_Comparison,
                plot_LDheatmap
            ))
        ) # visualization End
        #________________________________________

        import_p_link <- tabPanel(
            "p_link",
            tabPanelBody(
                "IMPORT_p_link",
                fluidRow(column(width = 3), column(width = 6, p("File"))),
                fluidRow(column(width = 3,
                                actionButton("ped_bt", label = "Choose ped File")),
                         column(width = 6,
                                fluidRow(textInput("ped_in",NULL)),
                                fluidRow(radioButtons("ped_sep", label = "Seprate of ped file",
                                                      choiceNames = c("tab","space"),
                                                      choiceValues = c("\t", " ")))),
                         column(width = 4,textOutput("ped_file"))),
                fluidRow(column(width = 3,
                                actionButton("map_bt", label = "Choose map File")),
                         column(width = 6,
                                fluidRow(textInput("map_in",NULL)),
                                fluidRow(radioButtons("map_sep", label = "Seprate of map file",
                                                      choiceNames = c("tab","space"),
                                                      choiceValues = c("\t", " ")))),
                         column(width = 4,textOutput("map_file"))),
                fluidRow(column(width = 3,
                                actionButton("plink_im", label = "Import"))),
                fluidRow(column(width = 12,
                                verbatimTextOutput("plink_info")))))

        import_fasta <- tabPanel(
            "Fasta_file",
            tabPanelBody(
                "IMPORT_FASTA",
                br(),
                fluidRow(column(width = 12,
                                textInput("seq_in", "File", width = "100%"))),
                fluidRow(column(width = 6,
                                actionButton("seq_bt",width = "100%", label = "Choose Fasta File")),
                         column(width = 6,
                                actionButton("seq_im",width = "95%", label = "Import"))),
                fluidRow(column(width = 12,
                                tableOutput("seq_data")))))
        import_vcf <- tabPanel(
            "VCF file",
            tabPanelBody(
                "IMPORT_VCF",
                br(),
                fluidRow(column(width = 12,
                                textInput("vcf_in","File", width = "100%"))),
                fluidRow(column(width = 6,
                                actionButton("vcf_bt",width = "100%", label = "Choose VCF File")),
                         column(width =6,
                                fluidRow(actionButton("vcf_im",width = "95%", label = "Import")))
                ),
                fluidRow(column(width = 12,
                                verbatimTextOutput("vcf_info")))))

        import_genomic <- tabPanel( # import Genotype start
            "Genotypic Data",

            navlistPanel(
                # import VCF
                import_vcf,

                # import fasta
                import_fasta,

                # import p_link
                import_p_link,
            )
        )

        import_annotation <- tabPanel(
            "Annotation",
            tabPanelBody(
                "import annotations",
                br(),
                fluidRow(column(width = 3,
                                actionButton("gff_bt", width = "100%",
                                             label = "Choose GFF/BED file"),
                                actionButton("gff_im", width = "100%",
                                             label = "Import")),
                         column(width = 5,
                                textInput("gff_in", "File")),
                         column(width = 3,
                                radioButtons("gff_format", label = "Format",
                                             choiceNames = c("GFF","BED"),
                                             choiceValues = c("GFF","BED")))),
                fluidRow(column(width = 12,
                                verbatimTextOutput("annotation_data"))))
        )

        import_pheno <- tabPanel(
            "Phenotypic Data",
            tabPanelBody(
                "Import Phenotypic Data",
                br(),
                fluidRow(column(width = 4,
                                actionButton("pheno_bt", width = "100%", label = "Choose phenotype file (txt/csv)"),
                                actionButton("pheno_im", width = "100%", label = "Import")),
                         column(width = 5,
                                textInput("pheno_in", label = "File")),
                         column(width = 3,
                                radioButtons("pheno_format",label = "format",
                                             choiceNames = c("txt","csv"),
                                             choiceValues = c("txt","csv")))),
                fluidRow(column(width = 12,
                                checkboxInput("pheno_showdata",
                                              label = "Show pheno data after import",
                                              value = FALSE))),
                fluidRow(column(width = 12,
                                htmlOutput("pheno_info"))),
                fluidRow(column(width = 12,
                                tableOutput("pheno_data"))))
        )

        import_infos <- tabPanel(
            "Accession Information Data",
            tabPanelBody(
                "import accession information data",
                br(),
                column(width = 6,
                       fluidRow(
                           textInput("accinfo_in", label = "File"),
                           radioButtons("accinfo_format",label = "Format",
                                        choiceNames = c("txt","csv"),
                                        choiceValues = c("txt","csv")),
                           actionButton("accinfo_bt", label = "Choose Accession information file"),
                           actionButton("accinfo_im", label = "Import"),
                           checkboxInput("accinfo_showdata",
                                         label = "Show accession information data after import",
                                         value = FALSE),
                           htmlOutput("accinfo_info"),
                           tableOutput("accinfo_data"))),
                column(width = 6,
                       fluidRow(
                           textInput("accinfo_in1", label = "File"),
                           radioButtons("accinfo_format1",label = "Format",
                                        choiceNames = c("txt","csv"),
                                        choiceValues = c("txt","csv")),
                           actionButton("accinfo_bt1", label = "Choose Accession information file"),
                           actionButton("accinfo_im1", label = "Import"),
                           htmlOutput("accinfo_info1"),
                           tableOutput("accinfo_data1")))
            )
        )

        tabPanel_import <- tabPanel(
            "Import Data",
            tabsetPanel(# Genotype UI
                import_genomic,

                # Annotation UI
                import_annotation,

                # Phenotype UI
                import_pheno,

                # Accession info UI
                import_infos
            )
        )
        #_________________________________________________

        tabPanel_Haplotyping <- tabPanel( # Do haplotyping Start
            "Haplotype",
            fluidPage(
                fluidRow(
                    column(width = 2,
                           br(),
                           radioButtons("hapresult_source",
                                        label = "Genotypic Data Format",
                                        choiceNames = c("vcf", "p.link", "Fasta", "Table", "Hapmap"),
                                        choiceValues = c("vcf", "p.link", "Fasta", "Table", "Hapmap"),
                                        selected = "vcf")),
                    column(width = 5,
                           h3("Haplotype Identification"),
                           textInput("happrefix", "Prefix of haplotype names", value = "H"),
                           checkboxInput("hyb_remove", "Remove heterozygote individuals", value = TRUE),
                           checkboxInput("na_drop", "Remove genotype unknown individuals", value = TRUE),
                           fluidRow(column(width = 6,
                                           actionButton("hapresult_bt", "1. Identify Haplotype")),
                                    column(width = 6,
                                           actionButton("hapsummary_bt", "2. Summary")))),
                    column(width = 5,
                           h3("Save Result"),
                           textInput("resultdir", "Directory"),
                           textInput("resultname", "File Name"),
                           h4("Save as: "),
                           column(width = 6,
                                  radioButtons("saveresult", label = NULL,
                                               choices = list("hapResult" = "hapResult"))),
                           column(width = 6,
                                  radioButtons("saveresult", label = NULL,
                                               choices = list("hapSummary" = "hapSummary"))),
                           column(width = 8,
                                  # shinyFiles::shinyDirChoose("lll","jjj"),
                                  # shinyFiles::shinyDirButton("shinychooseresultdir",
                                  #                            "Choose Result Directory","xuanzejieguo"),
                                  actionButton("chooseresultdir", "Choose Result Directory")),
                           column(width = 4,
                                  actionButton("download_hapresult", "Save")),
                    ),
                    br(),
                    fluidRow(column(width = 6,
                                    h4("Haplotype result"),
                                    verbatimTextOutput("hapresult", )),
                             column(width = 6,
                                    h4("Summary of haplotype result"),
                                    verbatimTextOutput("hapsummary"))),
                ),
            )
        )


        fluidPage(titlePanel("geneHapR packages"), # title
                  # Main Steps
                  useShinyjs(),
                  navbarPage("geneHapR-GUI",
                             id = "geneHapR",

                             # Welcome UI
                             tabPanel_welcome,

                             # Import data UI
                             tabPanel_import,

                             # preprocess UI
                             tabPanel_preprocess,

                             # Haplotype Identification
                             tabPanel_Haplotyping,

                             # Visualization
                             tabPanel_Visualization))
    }



    is.surfix <- function(x, surfix) {
        tolower(x)
        return(any(endsWith(x, surfix)))

    }
    ffilter <- matrix(ncol = 2, byrow = TRUE,
                      c("fasta \u6587\u4EF6 (*.fa, *.fasta)","*.fa;*.fatsa",
                        "vcf \u6587\u4EF6 (*.vcf, *.vcf.gz)", "*.vcf;*.vcf.gz",
                        "p.link \u6587\u4EF6 (ped & map)","*.map;*.ped",
                        "hapmap \u6587\u4EF6 (*.hmp)"," *.hmp",
                        "CSV\u9017\u53F7\u5206\u9694\u7684\u6587\u672C\u6587\u4EF6 (*.csv)", "*.csv",
                        "\u5236\u8868\u7B26\u5206\u9694\u7684\u6587\u672C\u6587\u4EF6 (*.txt)","*.txt",
                        "hapResult ( *.hapResult, *.txt)", "*.hapResult;*.txt",
                        "hapSummary ( *.hapSummary, *.txt)", "*.hapSummary;*.txt",
                        "GFF \u6587\u4EF6 (*.gff, *.gff3)", "*.gff;*.gff3",
                        "BED \u6587\u4EF6 (*.bed, *.bed4, *.bed6)","*.bed;*.bed4;*.bed6",
                        # 所有文件
                        "\u6240\u6709\u6587\u4EF6 (*.*)","*.*"))
    row.names(ffilter) <- c("fa", "vcf", "plink","hmp",
                            "csv","txt","hapResult","hapSummary",
                            "gff","bed","all")
}

