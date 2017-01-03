library(shiny)
library(shinyIncubator)
library(stringr)
library(RMySQL)
library(knitr)
library(xtable)
library(log4r)
loggerServer <- create.logger()
logfile(loggerServer) <- 'logs/server.log'
level(loggerServer) <- 'ERROR'

options(stringsAsFactors=F)
host<-'localhost'
version<-'v15_20140605'

disableActionButton <- function(id,session) {
  session$sendCustomMessage(type="jsCode",
                            list(code=paste0("$('#",id,"').prop('disabled',true)")) )
}

enableActionButton <- function(id,session) {
  session$sendCustomMessage(type="jsCode",
                            list(code=paste0("$('#",id,"').prop('disabled',false)")))
}

sendAlertMessage <- function(mess,session) {
  session$sendCustomMessage(type="jsCode",
                            list(code=paste0("alert('",mess,"')")))
}

# initShinyVars <- function(session) {
#   session$sendCustomMessage(type="jsCode",list(code="Shiny.unbindAll();
# Shiny.shinyapp.$inputValues.exploreTablePage = undefined;
# Shiny.shinyapp.$inputValues.exploreTableVariation = undefined;
# Shiny.bindAll();")
#   )
# }

shinyServer(function(input, output, session) {
  values <- reactiveValues()
  values$batch <- NULL
  values$sessionid <- gsub(':','',paste(format(Sys.time(), "%d%b%Y_%X"),sep='_'))
  values$iname <- ''
  values$inputFileType <- 'simple'
  # can be 'simple' (sample_id gene variation), 
  #        'clinical' (Worksheet  InvID	DNA No	Gene	Nomenclature	Het-hom	Class	Exon)
  #     or 'VCF'- to be implemented
  values$sessionschanged <- 0
  values$filesRead <- NULL
  values$filesLoaded <- NULL
  values$lastExploreButton <- 0
  values$lastBatchButton <- 0
  values$lastParseButton <- 0
  values$lastLoadButton <- 0
  values$lastSaveButton <- 0
  values$lastReturnButton <- 0
  values$lastMarkButton <- 0
  values$lastAddNoteButton <- 0
  values$lastSaveNotesButton <- 0
  values$lastExcludeDeep <- TRUE
  values$lastExcludeExtra <- TRUE
  values$lastExcludeCommon <- FALSE
  values$lastExcludeIndels <- FALSE
  values$selectedRow <- NULL
  values$contentType <- 'front'
  values$exploreTable <- NULL
  values$selectedVar <- NULL
  values$selectedVarBaseline <- NULL
  values$selectedVarClinical <- NULL
  values$selectedVarClinicalCurr <- NULL
  values$exploreTableiDL <- 10
  values$exploreTableiDS <- 0
  values$showSaveNotes <- FALSE
  values$showNewNote <- FALSE
  values$newNoteCurr <- ''
  values$offsetLimit <- 10
  values$variantReport <- NULL
  values$sources <- c('ALAMUT','EVS','TGP','TGP_PHASE3','ICR','HGMD','IARC','DMUDB','BIC','LOVD',
                      'EASTON','LINDOR','HOUDAYER','GUIDUGLI','WALKER','RMH','MUTTASTER',
                      'POLYPHEN2','DROST','BOUWMAN','UMD','CADD', 'SUSPECT')
  disableActionButton('parseButton',session)
  disableActionButton('loadButton',session)
  fixempty <- function(x) {return(ifelse(is.null(x)||is.na(x)||x==0,'',x))}
  freqprint <- function(x){return(sprintf('%5.4f',as.numeric(x)))}
  mydb<-dbConnect(MySQL(),user='anonymous',dbname='cigma2',host=host)
  rs<-dbSendQuery(mydb,"select table_name,field_name,field_label,variant_report,variant_report_order,analysis_output,analysis_output_order from fields where online_report='Y'")
  f<-fetch(rs,-1)
  fn<-split(f$field_name,f$table_name)
  fl<-split(f$field_label,f$table_name)
  # fa is a structure with the fields for the analysis output
  # the fields from the main classification and dogma_batchlog tables have to be entered manually
  fa<-f[which(f$analysis_output=='Y' & !(f$table_name %in% c('main','classification','dogma_batchlog'))),]
  rownames(fa)<-paste(gsub('tgp_phase3_pop','tgp_phase3',gsub('tgp_pop','tgp',fa[,1])),fa[,2],sep='_')
  fv<-f[which(f$variant_report=='Y' & !(f$table_name %in% c('classification','dogma_batchlog'))),]
  rownames(fv)<-paste(gsub('tgp_phase3_pop','tgp_phase3',gsub('tgp_pop','tgp',fv[,1])),fv[,2],sep='_')
  fan<-split(fa$field_name,fa$table_name)
  fvn<-split(fv$field_name,fv$table_name)
  
  # fields is a structure for the display panels
  # create list of tables with named lists of fields descriptions:
  fields<-setNames(lapply(names(fl),function(x){return(setNames(paste0(fl[[x]],':'),fn[[x]]))}),names(fl))
  # to access the description fieldX from tableY, use: fields[['tableX']][['fieldX']]
  
  #hack for tgp_pop & tgp_phase3_pop - these secondary tables store the information in multiple rows, one for each population
  fields[["tgp_pop"]]<-gsub(' .EUR.','',fields[["tgp_pop"]][1:4])
  fields[["tgp_pop"]]<-setNames(fields[["tgp_pop"]],gsub('_EUR','',names(fields[["tgp_pop"]])))
  fields[["tgp_phase3_pop"]]<-gsub(' .ACB.','',fields[["tgp_phase3_pop"]][1:4])
  fields[["tgp_phase3_pop"]]<-setNames(fields[["tgp_phase3_pop"]],gsub('_ACB','',names(fields[["tgp_phase3_pop"]])))
  
  rs<-dbSendQuery(mydb,"select distinct gene,ensembl transcript, substring_index(refseq,'.',1) rtranscript from cappagenes")
  #                select distinct gene,transcript,substring_index(rtranscript,'.',1) rtranscript from main")
  genelist<-fetch(rs,-1)
  gl<-as.vector(apply(genelist,1,function(x){paste(x[1]," (",x[2],")",sep="",collapse="")}))
  
  rs<-dbSendQuery(mydb,"select * from tgp_populations")
  dataset<-fetch(rs,-1)
  tp<-as.list(apply(dataset,1,function(x){ifelse(x[3]=='',x[2],paste(x[2]," (Super-population=",x[3],")",sep="",collapse=""))}))
  tp<-setNames(tp,dataset[,1])
  
  dbDisconnect(mydb)
  
  observe ({
    if (!is.null(input$offsetLimit) && input$offsetLimit>0)
      values$offsetLimit<-input$offsetLimit
  })
  
  # explore stuff
  currentgenetrans <- '' 
  # Initialize variationText every time we change the gene 
  observe({
    if (!is.null(input$genetrans) && (input$genetrans!=currentgenetrans || input$genetrans=='')) {
      updateTextInput(session, 'variationText', value='')
      currentgenetrans <- input$genetrans
    }
  })
  
  output$varCountText<-renderText({
    paste(varCount(),'variant(s)')
  })
  
  # Update the variation select list based on the selected gene 
  # and the initial input in the variation text box
  varCount <- reactive({
    splitgene<-unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))
    gene<-splitgene[1]
    transcript<-splitgene[2]
    variation<-input$variationText
    shortp <- substr(variation,4,4)>='0' && substr(variation,4,4)<='9'
    v<-NULL
    if (nchar(gene)>0 && nchar(variation)>3 && grepl('^[cp][.]',variation,perl=T)) {
      # don't bother unless there is a gene selected and at least 4 characters in variation
      mydb<-dbConnect(MySQL(),user='anonymous',dbname='cigma2',host=host)
      if (substr(variation,1,2) == 'c.') {
        rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                   "select CONCAT(hgvs_cdna,' (',hgvs_prot,')') from main_stable m where m.gene='",gene,
                                   "' and m.transcript='",transcript,"' and m.hgvs_cdna like '",variation,"%'",
                                   " order by cdna_pos,offset"))
      }else{ # it has to be 'p.'
        if (shortp) {# short p. notation
          rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                     "select CONCAT(hgvs_cdna,' (',hgvs_prot_code1,')') from main_stable m where m.gene='",gene,
                                     "' and m.transcript='",transcript,"' and m.hgvs_prot_code1 like '",variation,"%'",
                                     " order by cdna_pos,offset"))
        }else{ # long p notation
          rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                     "select CONCAT(hgvs_cdna,' (',hgvs_prot,')') from main_stable m where m.gene='",gene,
                                     "' and m.transcript='",transcript,"' and m.hgvs_prot like '",variation,"%'",
                                     " order by cdna_pos,offset"))
        }
      }
      if(!(grepl('^p[.]',variation,perl=T) && !shortp && nchar(variation)<6)) { 
        # don't do long p queries unless at least 6 char long
        v<-fetch(rs,-1)
        dbDisconnect(mydb)
        # Change values for input$variationSel if there are any results
        if(!is.null(v) && class(v)=='data.frame' && dim(v)[1]){
          updateSelectInput(session, "variationSel", choices = as.vector(v[,1]))
          return(dim(v)[1])
        }else{
          updateSelectInput(session, "variationSel", choices = c(''))
          return(0)
        }
      }else{
        v<-fetch(rs,-1)
        dbDisconnect(mydb)
        updateSelectInput(session, "variationSel", choices = c(''))
        return(0)
      }
    }else{
      updateSelectInput(session, "variationSel", choices = c(''))
      return(0)
    }
  })
  
  # Reactive function collecting loaded sets
  sesslist <- reactive({
    values$sessionschanged
    mydb<-dbConnect(MySQL(),user='anonymous',dbname='cigma2',host=host)
    rs<-dbSendQuery(mydb,"select investigator_name,session_id,count(*) variations from dogma_batch group by investigator_name,session_id")
    sl<-fetch(rs,-1)
    dbDisconnect(mydb)
    return(as.vector(apply(sl,1,function(x){paste(x[1]," (",x[2],", ",x[3]," variants)",sep="",collapse="")})))
  })
  
  # Dynamic UI with the loaded sets
  output$loadedSets<-renderUI({
    sl<-sesslist()
    info(loggerServer,paste("Loaded sets='",sl,"'"))
    if (values$sessionschanged == 0) {
      cs<-''
    }else{
      cs<-sl[grep(isolate(values$sessionid),sl)]
    }
    info(loggerServer,paste("Selected set='",cs,"'",sep=""))
    tagList(tags$script(src = "select2-master/select2.js"), 
            tags$link(href="select2-master/select2.css",rel="stylesheet"),
            tags$script(src = "js/include_select2.js"), 
            do.call(selectInput,list(inputId = "batchid", label="Previously loaded set :", choices = c('',sl), selected=cs, selectize=FALSE))
    )
  })
  
  ###########
  # BUTTONS #
  ###########
  
  # Test header of clinical lab file
  fileIsClinical <- function(x) {
    if (grepl('worksheet',x[1],ignore.case=T,perl=T) &&
          grepl('invid',x[2],ignore.case=T,perl=T) &&
          grepl('dna.no',x[3],ignore.case=T,perl=T) &&
          grepl('gene',x[4],ignore.case=T,perl=T) &&
          grepl('nomenclature',x[5],ignore.case=T,perl=T) &&
          grepl('het.hom',x[6],ignore.case=T,perl=T) &&
          grepl('class',x[7],ignore.case=T,perl=T) &&
          grepl('exon',x[8],ignore.case=T,perl=T)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  # Upload reactive function
  observe({
    # input$uploadFile will be NA initially. After the user selects and uploads a 
    # file, it will be a data frame with 'name', 'size', 'type', and 'datapath' 
    # columns. The 'datapath' column will contain the local filenames where the 
    # data can be found.
    inFile <- input$uploadFile
    error <- F
    if (!is.null(inFile) && !is.na(inFile) && (is.null(values$fileRead) || is.na(values$fileRead[inFile$datapath])))  {
      info(loggerServer,paste("Reading file:'",inFile$datapath, "'",sep=""))
      values$inputFileType <- 'simple'
      firstline<-scan(inFile$datapath, what=character(), sep="\t", quote='', nlines = 1)
      if (length(which(grepl('^gene$',firstline,ignore.case=T,perl=T)))) { # tab-delimited, no quotes
        if (fileIsClinical(firstline)){values$inputFileType <- 'clinical'}
        batch <- read.csv(inFile$datapath, header=TRUE, sep='\t', quote='')   #header=input$header, sep=input$sep, quote=input$quote)
      }else{
        firstline<-scan(inFile$datapath, what=character(), sep=",", quote='"', nlines = 1)
        if (length(which(grepl('^gene$',firstline,ignore.case=T,perl=T)))) { # comma-delimited, quotes='"'
          if (fileIsClinical(firstline)){values$inputFileType <- 'clinical'}
          batch <- read.csv(inFile$datapath, header=TRUE, sep=',', quote='"')   #header=input$header, sep=input$sep, quote=input$quote)
        }else{
          error(loggerServer,"Error - invalid input file")
          sendAlertMessage("Invalid input file! \\n\\nPlease use a tab-delimited file with no quotes\\nor a comma-delimited file with quotes (like the Excel .csv files)\\ncontaining at least these 3 fields specified on the first line:\\n\\nSampleID, Gene, Variant",session)
          error = T
        }
      }
      if (!error) {
        debug(loggerServer,paste("batch has these columns:",paste(colnames(batch),collapse=","),sep=""))
        debug(loggerServer,paste("batch has",dim(batch)[1],"rows"))
        values$fileRead[inFile$datapath]<-1
        values$batch<-batch
        values$contentType<-'table'
        values$sessionid <- gsub(':','',paste(format(Sys.time(), "%d%b%Y_%X"),sep='_'))
        enableActionButton('parseButton',session)
        disableActionButton('loadButton',session)
      }
    }
  })
  # Parse reactive function
  observe({
    if (is.null(input$parseButton) || is.na(input$parseButton) || input$parseButton == 0)
      return()
    if (input$parseButton != values$lastParseButton) {
      values$lastParseButton <- input$parseButton
      disableActionButton('parseButton',session)
      batch <- isolate(values$batch)
      if(!is.null(batch)) {
        progress <- Progress$new(session, min=1, max=10)
        on.exit(progress$close())
        progress$set(message = 'Parsing in progress',
                     detail = '  looking for sample_id, gene and variation')
        if (length(which(is.na(batch[,1])))) {batch <- batch[-which(is.na(batch[,1])),]}
        colnames(batch)[which(grepl('dna.no|sample',colnames(batch),ignore.case=T,perl=T))]='sample_id'
        colnames(batch)[which(grepl('gene',colnames(batch),ignore.case=T,perl=T))]='gene'
        colnames(batch)[which(grepl('nomenclature|variation|variant|hgvs',colnames(batch),ignore.case=T,perl=T))]='variation'
        debug(loggerServer,paste("After parsing batch has these columns:",paste(colnames(batch),collapse=","), sep=""))
        progress$set(detail = ' cleaning up the variation field', value=5)
        batch[,'variation']=sub(';.*$','',batch[,'variation'],perl=T)
        batch[,'variation']=sub('[_\\(]p.*$','',batch[,'variation'],perl=T)
        batch[,'variation']=sub(',c.*$','',batch[,'variation'],perl=T)
        values$batch<-batch#[which(grepl('BRCA',batch$gene)),]
        values$contentType<-'table'
        info(LoggerServer," Parsing done")
        progress$set(detail = ' done!', value=10)
        enableActionButton('loadButton',session)
      }
    }
  })
  
  
  # function generating SQL code for insert in database of a loaded set
  make_insert <- function(x) {
    i <- "INSERT INTO dogma_batch (session_id,investigator_name,sample_id,gene,variation"
    if (isolate(values$inputFileType == 'clinical')) {
      i <- paste0(i,',worksheet,inv_id,zygosity,exon')
    }
    i <- paste0(i,") VALUES ('",
                paste(isolate(values$sessionid),isolate(values$iname),x["sample_id"],x["gene"],x["variation"],sep="','"),
                "'"
    )
    if(isolate(values$inputFileType=='clinical')){
      i <- paste0(i,",'",
                  paste(x["Worksheet"],x["InvID"],x["Het.hom"],x["Exon"],sep="','"),
                  "'"
      )
    }
    i <- paste0(i,")")
    i
  }
  
  # Load & Analyze reactive function
  observe({
    if (is.null(input$loadButton) || is.na(input$loadButton) || input$loadButton == 0)
      return()
    if (input$loadButton != values$lastLoadButton) {
      if (isolate(input$iname=='')) {
        session$sendCustomMessage(type="jsCode",
                                  list(code="alert('Please fill in the investigator name!');$('input#iname').focus();"));
      }else{
        values$iname<-isolate(input$iname)
        values$lastLoadButton <- input$loadButton
        disableActionButton('loadButton',session)
        batch <- isolate(values$batch)
        if(!is.null(batch) && (is.null(values$fileLoaded) || is.na(values$fileLoaded[values$sessionid]))) {
          info(loggerServer," Inserting rows in the database")
          if (length(which(grepl('^$|no mutation',batch[,'variation'],ignore.case=T,perl=T)))) {
            batch <- batch[-which(grepl('^$|no mutation',batch[,'variation'],ignore.case=T,perl=T)),]
          }
          progress <- Progress$new(session, min=1, max=dim(batch)[1])
          on.exit(progress$close())
          progress$set(message = 'Loading in progress',
                       detail = paste('  Inserting',max=dim(batch)[1],'rows in the database'))
          mydb<-dbConnect(MySQL(),user='batch',password='cigma',dbname='cigma2',host=host)
          c <- 0
          for (i in apply(batch,1,make_insert)) {
            c <- c + 1
            progress$set(value = c)
            rs<-dbSendQuery(mydb,i)
          }
          rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                     "update dogma_batch b join cappagenes c 
             on b.gene=c.gene 
            set b.transcript=c.ensembl 
          where b.transcript=''
            and b.session_id ='",isolate(values$sessionid),"'"))
          info(loggerServer," Analyzing variants in the database")
          procs<-c("CALL get_alleles_from_hgvs_cdna('dogma_batch','variation')",
                   "CALL get_varType_from_hgvs_cdna('dogma_batch','variation')",
                   "CALL get_varLocation_from_hgvs_cdna('dogma_batch','variation')",
                   "CALL get_cdnacoord_from_hgvs_cdna('dogma_batch','variation')",
                   "CALL get_genomecoord_from_cdnacoord('dogma_batch')",
                   "CALL determine_codingEffect('dogma_batch')",
                   "CALL make_hgvs('dogma_batch')")
          for (i in procs) {
            progress$set(message = 'Analyzing in progress',detail = unlist(str_split(i,'[(]'))[1])
            rs<-dbSendQuery(mydb,i)
          }
          info(loggerServer," Analyzing done")
          dbDisconnect(mydb)
          values$fileLoaded[values$sessionid]<-1      
          values$sessionschanged<-values$sessionschanged + 1
        }
      }
    }
  })
  
  # Select an already loaded set reactive function
  observe({
    x<-input$batchid
    debug(loggerServer,paste("Select loaded set changed to '",x,"'",sep=""))
    if (!is.null(x) && !is.na(x)) {
      if (x !='') {
        progress <- Progress$new(session, min=1, max=10)
        on.exit(progress$close())
        progress$set(message = 'Loading in progress',
                     detail = ' querying database')
        sessid<-unlist(str_split(unlist(str_split(x,' [(]'))[2],'[,]'))[1]
        iname<-unlist(str_split(x,' [(]'))[1]
        values$iname <- iname
        values$sessionid <- sessid
        debug(loggerServer,paste(" sessid=",sessid," iname=",iname))
        exploreTable<-NULL
        mydb<-dbConnect(MySQL(),user='batch',password='cigma',dbname='cigma2',host=host)
        rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                   "create temporary table tmp_",sessid," as 
                                   select gene,variation,
                                          max(concat_ws('\t',date_format(creat_date,'%Y-%m-%d'),clinical_cigma_class)) last_clinical,
                                          group_concat(concat_ws('\t',investigator_name,clinical_cigma_class,
                                                                 concat(date_format(creat_date,'%d/%m/%Y'),' (',
                                                                        if(round(datediff(curdate(),creat_date)/7)<9,
                                                                           concat(round(datediff(curdate(),creat_date)/7),' weeks ago'),
                                                                           concat(round(datediff(curdate(),creat_date)/30.4166),' months ago')
                                                                          ),')' ),
                                                                  notes
                                                                 )
                                                         order by creat_date desc separator '\n'
                                                        ) notescontent 
                             from dogma_batchlog group by gene,variation"));
        rs<-dbSendQuery(mydb,paste(sep="",collapse="","create index i1 on tmp_",sessid," (gene,variation)"))
        rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                   "select b.sample_id,
                IF(b.hgvs_prot!='',CONCAT(b.gene,':',b.variation,' (',IF(b.codingEffect='synonymous','p.=',b.hgvs_prot),')')
                                  ,CONCAT(b.gene,':',b.variation)) variation,
                c.cigma_class working_cigma_class,
                c.classification_justification,
                IF(b.clinical_cigma_class!='', b.clinical_cigma_class, substring_index(d.last_clinical,'\t',-1)) clinical_cigma_class,
                IF(b.notes!='',substring_index(substring_index(b.notes,'\t',3),'\t',-1),
                               date_format(str_to_date(substring_index(d.last_clinical,' ',1),'%Y-%m-%d'),'%d/%m/%Y')) clinical_cigma_date,
                d.notescontent, b.notes notes_new, 
                b.varLocation, b.varType, b.codingEffect, b.offset,
                b.worksheet, b.inv_id, b.zygosity, b.exon, case when cv.reason is not null then 1 else 0 end common
           from dogma_batch b left join dogma_classification c
             on b.gene=c.gene and b.variation=c.hgvs_cdna and c.version='",version,"'
                left join tmp_",sessid," d
             on b.gene=d.gene and b.variation=d.variation
                left join commonvar cv
             on b.gene=cv.gene and b.variation=cv.hgvs_cdna
           where session_id ='",sessid,"'"))
        exploreTable<-fetch(rs,-1)
        exploreTable[is.na(exploreTable[,3]),3]<-''
        exploreTable[is.na(exploreTable[,4]),4]<-''
        exploreTable[is.na(exploreTable[,5]),5]<-''
        exploreTable[is.na(exploreTable[,6]),6]<-''
        exploreTable[is.na(exploreTable[,7]),7]<-''
        exploreTable[is.na(exploreTable[,8]),8]<-''
        exploreTable[is.na(exploreTable[,13]),13]<-''
        exploreTable[is.na(exploreTable[,14]),14]<-''
        exploreTable[is.na(exploreTable[,15]),15]<-''
        exploreTable[is.na(exploreTable[,16]),16]<-''
        dbDisconnect(mydb)
        progress$set(detail =' creating table', value  = 5)
        exploreTable<-cbind(notes='', exploreTable)
        exploreTable[union(which(exploreTable$notescontent!=''),which(exploreTable$notes_new!='')),'notes']<-'<img src="images/details_open.png">'
        exploreTable$visited <-''
        debug(loggerServer,paste(" exploreTable has: ",dim(exploreTable)," dimensions"))
        values$exploreTable<-exploreTable
        values$inputFileType<-ifelse(exploreTable[1,14]=='','simple','clinical')
        values$contentType<-'exploreTable'
        disableActionButton('parseButton',session)
        disableActionButton('loadButton',session)
        progress$set(detail =' done!', value = 10)
      }else{
        values$batch<-NULL
        values$contentType<-'table'
        disableActionButton('parseButton',session)
        disableActionButton('loadButton',session)
      }
    }
  })
  
  
  # function generating SQL code for update of the notes field
  make_update_notes <- function(x) {
    paste0("UPDATE dogma_batch SET notes='",x["notes_new"],"' WHERE session_id='",x["session_id"],"' AND variation='",x["variation"],"'" )
  }
  
  # function generating SQL code for update of the clinical_cigma_class field
  make_update_clinical <- function(x) {
    paste0("UPDATE dogma_batch SET clinical_cigma_class='",x["clinical_cigma_class"],"' WHERE session_id='",x["session_id"],"' AND variation='",x["variation"],"'")
  }
  
  # saveButton reactive function
  observe({
    info(loggerServer,"in function for saveButton")
    if (is.null(input$saveButton) || is.na(input$saveButton) || input$saveButton == 0)
      return()
    debug(loggerServer,paste("Save button pressed:",input$saveButton,sep=""))
    exploreTable <- isolate(values$exploreTable)
    if(!is.null(exploreTable)) {
      x<-isolate(input$batchid)
      debug(loggerServer,paste(" Saving batch",x))
      if (length(which(exploreTable$clinical_cigma_class!=''))+length(which(exploreTable$notes_new!=''))>0) {
        batch<-unique(exploreTable[union(which(exploreTable$clinical_cigma_class!=''),which(exploreTable$notes_new!='')),c('variation','clinical_cigma_class','notes_new')])
        batch$session_id<-unlist(str_split(unlist(str_split(x,' [(]'))[2],'[,]'))[1]
        batch$variation<-sub(' .*$','',sub('^.*:','',batch$variation,perl=T),perl=T)
        rows1<-length(which(batch$notes_new!=''))
        rows2<-length(which(batch$clinical_cigma_class!=''))
        progress <- Progress$new(session, min=1, max=rows1+rows2)
        on.exit(progress$close())
        progress$set(message = 'Saving in progress',
                     detail = paste('  Updating',max=rows1+rows2,'rows in the database'))
        mydb<-dbConnect(MySQL(),user='batch',password='cigma',dbname='cigma2',host=host)
        c <- 0
        for (i in apply(batch[which(batch$notes_new!=''),],1,make_update_notes)) {
          c <- c + 1
          progress$set(detail = paste('  Updating notes for:',batch$variation[c]),value = c)
          rs<-dbSendQuery(mydb,i)
        }
        for (i in apply(batch[which(batch$clinical_cigma_class!=''),],1,make_update_clinical)) {
          c <- c + 1
          progress$set(detail = paste('  Updating curated class for:',batch$variation[c-rows1]),value = c)
          rs<-dbSendQuery(mydb,i)
        }
      }
      info(loggerServer," Saving done")
      dbDisconnect(mydb)
      values$batch<-NULL
      values$sessionschanged<-0
      values$lastExcludeDeep <- TRUE
      values$lastExcludeExtra <- TRUE
      values$lastExcludeCommon <- FALSE
      values$lastExcludeIndels <- FALSE
      values$selectedVar <- NULL
      values$selectedVarBaseline <- NULL
      values$selectedVarClinical <- NULL
      values$batchid<-NULL
      values$exploreTable<-NULL
      values$exploreTableiDL<-10
      values$exploreTableiDS<-0
      values$contentType<-'table'
    }
  })
  
  # Explore table click resulting in input$exploreTableVariation change
  observe({
    if(!is.null(input$exploreTableVariation) && input$exploreTableVariation!='') {
      info(loggerServer,"in function for click on variation.")
      x <- unlist(str_split(input$exploreTableVariation, "\t"))
      if (length(x)){
        values$selectedVar <- gsub(' .*$','',x[1],perl=T)
        values$selectedVarBaseline <- x[2]
        values$selectedVarClinical <- x[3]
        values$selectedVarClinicalCurr <- x[3]
        values$exploreTableiDS <- as.numeric(x[4])
      }
      info(loggerServer,paste(" Selected variation:",isolate(values$selectedVar),sep=""))
      isolate({
        if(length(which(values$exploreTable$visited!=''))) {
          values$exploreTable[which(values$exploreTable$visited!=''),'visited']<-'Y'
          debug(loggerServer,paste("these rows have been set to Y:",which(values$exploreTable$visited!='')))
        }
        values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'visited']<-'L'
        debug(loggerServer,paste("these rows have been set to L:",which(grepl(values$selectedVar,values$exploreTable$variation))))
        values$newNoteCurr<-unlist(str_split(values$exploreTable[which(values$exploreTable$variation==x[1]),'notes_new'][1],"\t"))[4]
        if(is.null(values$newNoteCurr) || is.na(values$newNoteCurr)){values$newNoteCurr=''}
      })
      values$contentType <- 'browseTabs'
    }
  })
  
  # Explore table change length resulting in input$exploreTableiDSiDL change
  observe({
    if(!is.null(input$exploreTableiDSiDL) && input$exploreTableiDSiDL!='') {
      info(loggerServer,"in function for change info on exploreTable")
      x <- unlist(str_split(input$exploreTableiDSiDL, "\t"))
      if (length(x)){
        values$exploreTableiDS <- as.numeric(x[1])
        values$exploreTableiDL <- as.numeric(x[2])
      }
    }
  })
  
  # Explore (undo icon) reactive function
  observe({
    debug(loggerServer,paste("in function for exploreButton button=",input$exploreButton," last=",isolate(values$lastExploreButton)))
    if (is.null(input$exploreButton) || is.na(input$exploreButton)) {
      return()
    }
    if (input$exploreButton == 0) {
      values$lastExploreButton<-0
      return()
    }
    if ((input$exploreButton != isolate(values$lastExploreButton))) {
      values$lastExploreButton <- input$exploreButton
      values$contentType <- 'exploreSingle'
    }
  })
  
  # Batch (undo icon) reactive function
  observe({
    debug(loggerServer,paste("in function for batchButton button=",input$batchButton," last=",isolate(values$lastBatchButton)))
    if (is.null(input$batchButton) || is.na(input$batchButton)) {
      return()
    }
    if (input$batchButton == 0) {
      values$lastBatchButton<-0
      return()
    }
    if ((input$batchButton != isolate(values$lastBatchButton))) {
      values$lastBatchButton <- input$batchButton
      values$contentType <- 'table'
    }
  })
  
  # Return (undo icon) reactive function
  observe({
    debug(loggerServer,paste("in function for returnButton button=",input$returnButton," last=",isolate(values$lastReturnButton)))
    if (is.null(input$returnButton) || is.na(input$returnButton)) {
      return()
    }
    if (input$returnButton == 0) {
      values$lastReturnButton<-0
      return()
    }
    if ((input$returnButton != isolate(values$lastReturnButton))) {
      values$lastReturnButton <- input$returnButton
      if(values$contentType=='browseTabs'){
        if (isolate(values$showSaveNotes)){saveNotes()}
        values$newNoteCurr <- ''
        values$showNewNote <- FALSE
        values$showSaveNotes <- FALSE
        values$contentType <- 'exploreTable'
      }else{
        values$contentType <-'front'
      }
    }
  })
  
  # Mark variation reactive function
  observe({
    debug(loggerServer,paste("in function for markButton button=",input$markButton," last=",isolate(values$lastmarkButton)))
    if (is.null(input$markButton) || is.na(input$markButton)) {
      return()
    }
    if (input$markButton == 0) {
      values$lastmarkButton<-0
      return()
    }
    if ((input$markButton != isolate(values$lastmarkButton))) {
      values$lastmarkButton <- input$markButton
      if (input$clinicalCIGMA==''){
        sendAlertMessage("Please select a curated class!",session)
      }else{
        values$newNoteCurr <- isolate({paste0(ifelse(values$newNoteCurr!='',values$newNoteCurr,''),
                                              'REVIEWED on ',format(Sys.time(), '%d/%m/%Y\n'))})
        values$showNewNote <- TRUE
        values$showSaveNotes <- TRUE
      }
    }
  })
  
  # Clinical CIGMA select change function
  observe({
    x<-input$clinicalCIGMA
    if (!is.null(x) && !is.na(x)) {
      info(loggerServer,paste("Curated class changed to '",x,"'",sep=""))
      if (x != isolate(values$selectedVarClinical)) {
        values$newNoteCurr <- isolate({paste0(ifelse(values$newNoteCurr!='',values$newNoteCurr,''),
                                              'Curated class changed to ',x,'\n')})
        values$showNewNote <- TRUE
        values$showSaveNotes <- TRUE
      }
      values$selectedVarClinicalCurr <- x
    }
  })
  
  # Add Note reactive function
  observe({
    debug(loggerServer,paste("in function for addNoteButton button=",input$addNoteButton," last=",isolate(values$lastaddNoteButton)))
    if (is.null(input$addNoteButton) || is.na(input$addNoteButton)) {
      return()
    }
    if (input$addNoteButton == 0) {
      values$lastaddNoteButton<-0
      return()
    }
    if ((input$addNoteButton != isolate(values$lastaddNoteButton))) {
      values$showNewNote <- TRUE
      values$showSaveNotes <- TRUE
    }
  })
  
  saveNotes <- function() {
    isolate({
      values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'notes_new'] <- 
        paste(values$iname,
              input$clinicalCIGMA,
              format(Sys.time(), '%d/%m/%Y'),
              gsub('\n','<br>',input$newNote,perl=T),
              sep="\t"
        )
      values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'notes'] <- '<img src="images/details_open.png">'
      values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'clinical_cigma_class'] <- input$clinicalCIGMA
      values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'clinical_cigma_date'] <- format(Sys.time(), '%d/%m/%Y')
      values$newNoteCurr<-gsub('\n','<br>',input$newNote,perl=T)
    })
    values$showNewNote <- FALSE
    values$showSaveNotes <- FALSE
  }
  
  # Save Notes reactive function
  observe({
    debug(loggerServer,paste("in function for saveNotesButton button=",input$saveNotesButton," last=",isolate(values$lastSaveNotesButton)))
    if (is.null(input$saveNotesButton) || is.na(input$saveNotesButton)) {
      return()
    }
    if (input$saveNotesButton == 0) {
      values$lastsaveNotesButton<-0
      return()
    }
    if ((input$saveNotesButton != isolate(values$lastsaveNotesButton))) {
      saveNotes()
    }
  })
  
  
  ###################
  # DYNAMIC CONTENT #
  ###################
  # Output holding the dynamic content for the control panel
  output$uiControl<-renderUI({
    debug(loggerServer,"in renderUI for uiControl")
    c<-NULL
    if (values$contentType == 'front') {
      c<-div(class="container-fluid navbar navbar-fixed-top",
             div(class="well",style="margin-bottom:0;;padding:4px;",
                 fluidRow(div(align="center",img(src="images/CaVaDa.png"))
                 )
             )
      )
    }
    else if (values$contentType == 'exploreSingle') {
      c<-div(class="container-fluid navbar navbar-fixed-top",
             div(class="well",style="margin-bottom:0;;padding:4px;",
                 fluidRow(
                   column(width=2,img(src="images/dogma_logo_150.png")
                   ),
                   column(width=10,
                          fluidRow(
                            column(width=7, fluidRow(column(width=6,selectInput(inputId = "genetrans", label="Gene", choices = c('',gl), selectize=FALSE))),
                                   fluidRow(div(class="span6 thin","Variant starting with:"),
                                            div(class="span6 thin","Select matching variant:")),
                                   fluidRow(column(width=6,textInput(inputId = "variationText", label="")),
                                            column(width=6,selectInput(inputId = "variationSel", label="", choices = c(''), selectize=FALSE))),
                                   fluidRow(column(offset=6,width=6,textOutput("varCountText")))
                            ),
                            column(width=2, helpText(em('Examples: ')),
                                   helpText('BRCA1 c.181T>G'),
                                   #    helpText('BRCA1 c.557C>A'),
                                   helpText('BRCA2 c.7994A>G'),
                                   #    helpText('BRCA1   p.Arg1203Ter'),
                                   helpText('BRCA1 p.Y130X'),
                                   helpText('MLH1 c.131C>T'),
                                   helpText('MSH2 c.2500G>A')
                            ),
                            column(width=2,
                                   helpText(em("Classes:")),
                                   helpText("1-Non-pathogenic"),
                                   helpText("2-Likely non-pathogenic"),
                                   helpText("3-Uncertain/review"),
                                   helpText("4-Likely pathogenic"),
                                   helpText("5-Pathogenic")
                            ),
                            column(width=1,
                                   span(title="Return to main page",actionButton("returnButton","",icon=icon("home","fa-2x"))),
                                   span(title="Download PDF with variant report",downloadButton("downloadVariantReport",""))
                            )
                          )
                   )
                 )
             )
      )
    }
    else if (values$contentType == 'table') {
      c<-div(class="container-fluid navbar navbar-fixed-top",
             div(class="well",style="margin-bottom:0;;padding:4px;",
                 fluidRow(
                   column(width=2,img(src="images/dogma_logo_150.png")
                   ),
                   column(width=10,
                          fluidRow(
                            column(width=4, textInput(inputId = "iname", label="Investigator name:",value=values$iname)
                            ),
                            column(width=3, fileInput('uploadFile', 'Upload variants:',
                                                      accept=c('text/csv','text/comma-separated-value,text/plain','application/csv','.csv')),
                                   actionButton('parseButton',"Parse"),
                                   actionButton('loadButton',"Load & Analyze")
                            ),
                            column(width=4, uiOutput('loadedSets')
                            ),
                            column(width=1,
                                   span(title="Return to main page",actionButton("returnButton","",icon=icon("home","fa-2x"))
                                   )
                            )
                          )    
                   )
                 )
             )
      )
    }else if(values$contentType == 'exploreTable') {
      c<-div(class="container-fluid navbar navbar-fixed-top",
             div(class="well",style="margin-bottom:0;padding:4px;",
                 fluidRow(
                   column(width=2,img(src="images/dogma_logo_150.png")),
                   column(width=10,
                          fluidRow(
                            column(width=4,
                                   div(style="font-style:italic;font-size:17.5px;margin-bottom:10px","Set:"),
                                   div(style="font-weight:bold;font-size:17.5px;",textOutput("selectedSet")),
                                   div(style="margin-top:20px;",
                                       span(title="Return to main page",actionButton("returnButton","",icon=icon("home","fa-2x"))),
                                       span(title="Save changes made to the set",actionButton("saveButton","",icon=icon("save","fa-2x"))),
                                       span(title="Download the annotated set",downloadButton("downloadAnnotatedBatch",""))
                                   )
                            ),
                            column(width=5,
                                   checkboxInput("excludeDeep","Exclude intronic variants ",values$lastExcludeDeep),
                                   do.call(numericInput, list("offsetLimit","more than:",isolate(values$offsetLimit))),span(style="width:200px;","nt away from an exon border"),
                                   checkboxInput("excludeExtra","Exclude UTRs and intergenic",values$lastExcludeExtra),
                                   checkboxInput("excludeCommon","Exclude common variants",values$lastExcludeCommon),
                                   checkboxInput("excludeIndels","Exclude insertions/deletions",values$lastExcludeIndels)
                            ),
                            column(width=3,
                                   helpText("1-Non-pathogenic"),
                                   helpText("2-Likely non-pathogenic"),
                                   helpText("3-Uncertain/review"),
                                   helpText("4-Likely pathogenic"),
                                   helpText("5-Pathogenic")
                            )
                          )
                   )
                 )
             )
      )
    }else if (values$contentType == 'browseTabs') {
      saveNotes <- ''
      if(values$showSaveNotes){ 
        saveNotes<-span(title=paste("Save changes made to variant",isolate(values$selectedVar)),
                        do.call(actionButton,list(inputId="saveNotesButton",label="",icon=icon("save","fa-2x"))))
      }
      newNote <- span(title=paste("Add note for variant",isolate(values$selectedVar)),actionButton("addNoteButton","",icon=icon("plus-square-o","fa-2x")))
      if(values$showNewNote) {
        newNote <- tags$textarea(id="newNote", rows=6, cols=70, style="margin: 0px 0px 10px; width: 95%; height: 130px; resize: none", gsub('<br>','\n',values$newNoteCurr))
      }
      c<-div(class="container-fluid navbar navbar-fixed-top",
             div(class="well",style="margin-bottom:0;padding:4px;",
                 fluidRow(
                   column(width=2,img(src="images/dogma_logo_150.png")),
                   column(width=2,
                          div(style="font-style:italic;font-size:17.5px;","Variant:"),
                          div(style="margin-top:5px;font-weight:bold;font-size:17.5px;",textOutput("selectedVar")),
                          div(style="margin-top:20px;",
                              span(title=paste("Return to set:",isolate(input$batchid)),actionButton("returnButton","",icon=icon("reply","fa-2x"))),
                              span(title=paste("Mark variant",isolate(values$selectedVar),'as "REVIEWED on',format(Sys.time(), '%d/%m/%Y"')),actionButton("markButton","",icon=icon("check-square-o","fa-2x"))),
                              saveNotes
                          )
                   ),
                   column(width=2,
                          tags$style("#selectedVarBaseline {padding:4px;margin-bottom:5px;} #clinicalCIGMA {margin-bottom:5px;}"),
                          "Automated class:", verbatimTextOutput("selectedVarBaseline"),
                          selectInput("clinicalCIGMA","Curated class:",choices=c("",
                                                                                 "1-Non-pathogenic",
                                                                                 "2-Likely non-pathogenic",
                                                                                 "3-Uncertain",
                                                                                 "4-Likely pathogenic",
                                                                                 "5-Pathogenic"
                          ),
                          selected=values$selectedVarClinicalCurr, selectize=FALSE)
                   ),
                   column(width=3, newNote
                   ),
                   column(width=2,
                          helpText("1-Non-pathogenic"),
                          helpText("2-Likely non-pathogenic"),
                          helpText("3-Uncertain/review"),
                          helpText("4-Likely pathogenic"),
                          helpText("5-Pathogenic")
                   ),
                   column(width=1,
                          span(title="Download PDF with variant report",downloadButton("downloadVariantReport",""))
                   )
                 )
             )
      )
    }
    return(tagList(tags$script(src = "js/resize_download.js"),c)) 
  })
  
  # Output holding the dynamic content for the main panel
  output$uiContents<-renderUI({
    debug(loggerServer,"in renderUI for output$uiContent")
    debug(loggerServer,paste(" values$contentType='",values$contentType,"'",sep=''))
    if (values$contentType == 'front') {
      div(class="well",style="padding-top:170px;",
          fluidRow(
            column(width=3, img(src="images/dogma_logo_256.png")
            ),
            column(width=3, 
                   fluidRow(
                     div(align="center", style="margin-top:48px;",
                         helpText('Explore a single variant'),
                         actionButton('exploreButton',"Explore"),
                         div(style="margin-top:16px;",helpText('Explore a batch of variants'),
                             actionButton('batchButton',"Batch")
                         )
                     )
                   )
            ),
            column(width=6,
                   div(style="margin-top:6px;",
                       helpText(h4('Instructions:')),
                       helpText('- In ',strong('[Explore]'),' mode, select a gene and type the HGVS name of a single variant you wish to query (beginning with "c." or "p.")'),
                       helpText('- In ',strong('[Batch]'),' mode, upload a tab-delimited .txt or comma-separated .csv file with a header line including: ',em('Sample, Gene, Variant')),
                       fluidRow(column(offset=1, width=11, 
                                       helpText('- Click on ',strong('[Parse]'),' to process and ',strong('[Load&Analyze]'),' to load it in the database. A list of previously loaded sets can be viewed via the drop-down: any of these can be selected for re-examination.'),
                                       helpText('- Click on one variant name in the table for detailed information from related data sources'))),
                       helpText('- Click on ',img(src="images/details_open.png"),' for clinician curated notes on variant'),
                       helpText('- The ',strong('automated_class'),' is generated by the gene-specific decision-tree, with an associated ',strong('classification_justification'),'. The ',strong('curated_class'),' is assigned following clinician review +/- literature search.')
                   )
            )
          ),
          fluidRow(div(style="margin-top:6px;",
                       helpText(h4('Disclaimer:')),
                       helpText(em(p('The Cancer Variant Database is a repository of annotations for variants in cancer predisposition genes',br(),'Classifications, both automated and curated, are derived from integration of the variant level data according to objective criteria.'),
                                   p('This resource is currently under development and is for research use only.')))
          )
          )
      )
    }
    else if (values$contentType == 'table' && !is.null(values$batch)) {
      debug(loggerServer,paste(" values$batch columns:",paste(colnames(values$batch),collapse=","),sep=''))
      debug(loggerServer,paste(" values$batch has ",dim(values$batch)[1]," rows.",sep=''))
      tagList(
        singleton(
          tags$head(
            tags$style(type="text/css","div.dataTables_length label {float: right !important; text-align: right !important;}
                                    div.dataTables_info {float: right ! important;}
                                    div.dataTables_filter label {float: left !important;}
                                    div.dataTables_paginate {float: left !important; margin: 0;}"
            ))),
        div(class="well",style="padding-top:170px;",dataTableOutput('batch'))
      )
    }else if (values$contentType == 'exploreTable') {
      debug(loggerServer,paste(" exploreTable has:",dim(values$exploreTable)," dimensions"))
      div(class="well",style="padding-top:170px;",
          selDataTableOutput("exploreTable")
      )
    }else if (values$contentType == 'browseTabs' || values$contentType == 'exploreSingle') {
      debug(loggerServer,paste("selected variation=",values$selectedVar))
      info <- cigmainfo()
      notes_new <- NULL
      if  (values$contentType == 'browseTabs'){
        notes_new <- values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'notes_new']
      }
      tabs <- vector(mode="list")
      if (!is.null(info) && !is.null(info$summary)) {
        tabs$summary <- tabPanel("Summary", uiOutput('summaryPanel'))
      }
      if (!is.null(info) && 
            (!is.null(info$alamut) || 
               !is.null(info$muttaster) || 
               !is.null(info$polyphen2) || 
               !is.null(info$cadd) || 
               !is.null(info$suspect)
            )
      ) {
        tabs$insilico <- tabPanel("In silico predictions", uiOutput('insilicoPanel'))
      }
      if (!is.null(info) && 
            (!is.null(info$evs) || 
               !is.null(info$tgp) ||
               !is.null(info$tgp_phase3) ||
               !is.null(info$icr) ||
               !is.null(info$bic) ||
               !is.null(info$lovd) ||
               !is.null(info$dmudb) ||
               !is.null(info$umd)
            )
      ) {
        tabs$frequency <- tabPanel("Frequency", uiOutput('frequencyPanel'))
      }
      if (!is.null(info) && 
            (!is.null(info$easton) || 
               !is.null(info$lindor) ||
               !is.null(info$iarc) ||
               (!is.null(info$lovd) && info$lovd['insight_class',]!='')
            )
      ) {
        tabs$geneticepi <- tabPanel("Genetic Epidemiology", uiOutput('geneticepiPanel'))
      }
      if (!is.null(info) && 
            (!is.null(info$houdayer) || 
               !is.null(info$walker)
            )
      ) {
        tabs$splicing <- tabPanel("Splicing analysis", uiOutput('splicingPanel'))
      }
      if (!is.null(info) && 
            (!is.null(info$guidugli) || 
               !is.null(info$drost) || 
               !is.null(info$bouwman)
            )
      ) {
        tabs$functional <- tabPanel("Functional analysis", uiOutput('functionalPanel'))
      }
      notestabs <- vector(mode="list")
      if (!is.null(info)) { #&& !is.null(info$notes)) || (!is.null(notes_new) &&!is.null(notes_new[1]) && notes_new[1]!='')) {
        notestabs$notes <- tabPanel("CaVaDa classifications/notes", uiOutput('notesPanel'))
      }
      
      div(class="well",style="padding-top:170px;",
          fluidRow(
            column(width=8,do.call(tabsetPanel, c(id="tabpanel", tabs))),
            column(width=4,do.call(tabsetPanel, c(id="notespanel", notestabs)))
          )
      )
    }else{
      NULL
    }
  })
  
  #################################################
  # 1. Plain dataTable with loaded/parsed content #
  #    contentType = "table"                      #
  #################################################
  output$batch<-renderDataTable({
    values$batch
  }, options = function () {       
    list(sDom = "<'row-fluid'<'span6'p><'span6'i>>t<'row-fluid'<'span6'f><'span6'l>r>", 
         bDestroy = TRUE, 
         bLengthChange = TRUE, 
         iDisplayLength = values$exploreTableiDL,
         aLengthMenu= list(10, 25, 50, 100),
         fnInfoCallback = I("(function( oSettings, iStart, iEnd, iMax, iTotal, sPre ) {Shiny.onInputChange('exploreTableiDSiDL',oSettings._iDisplayStart+'\t'+oSettings._iDisplayLength);return(sPre);})")
    )}
  )
  
  ####################################
  # 2. dataTable with hidden fields  #
  #    contentType = "exploreTable"  #
  ####################################
  
  output$exploreTable <- renderDataTable({
    v<-values$exploreTable
    info(loggerServer,"in renderDataTable for output$exploreTable")
    debug(loggerServer,paste(" exploreTable has:",dim(v)," dimensions"))
    if (is.null(v)) return(v);
    if(input$excludeDeep && length(which(v$varLocation == 'intron' & abs(v$offset)>input$offsetLimit))) {
      v<-v[-which(v$varLocation == 'intron' & abs(v$offset)>input$offsetLimit),]
    }
    if(input$excludeExtra && length(which(v$varLocation %in% c("3'UTR","5'UTR","intergenic")))) {
      v<-v[-which(v$varLocation %in% c("3'UTR","5'UTR","intergenic")),]
    }
    if(input$excludeCommon && length(which(as.logical(v$common)))) {
      v<-v[-which(as.logical(v$common)),]
    }
    if(input$excludeIndels && length(which(v$varType %in% c('insertion','deletion','complex')))) {
      v<-v[-which(v$varType %in% c('insertion','deletion','complex')),]
    }
    colnames(v)[1]<-'Notes'
    colnames(v)[2]<-'SampleID'
    colnames(v)[3]<-'Gene:Variation'
    colnames(v)[4]<-'Automated class'
    colnames(v)[5]<-'Justification (automated class)'
    colnames(v)[6]<-'Curated class'
    colnames(v)[7]<-'Date'
    v
  },options = function () {
    list(sDom = "<'row-fluid'<'span5'p><'span3'T><'span4'i>>t<'row-fluid'<'span5'f><'span3'T><'span4'l>r>", 
         bDestroy = TRUE, 
         bLengthChange = TRUE, 
         iDisplayLength = values$exploreTableiDL,
         iDisplayStart = values$exploreTableiDS,
         aLengthMenu= list(10, 25, 50, 100),
         bSortClasses = TRUE,
         bAutoWidth = FALSE,
         aoColumnDefs = list(list(bVisible=FALSE, aTargets=sapply(c(7:18),list)),
                             list(sWidth='275px', aTargets=sapply(2,list)),
                             list(sWidth='50px' , aTargets=sapply(0:1,list)),
                             list(sWidth='225px' , aTargets=sapply(4,list)),
                             list(sWidth='160px', aTargets=sapply(c(3,5),list))),
         fnRowCallback  = I("(function (nRow,aData,iDisplayIndex,iDisplayIndexFull){if(aData[18]=='Y'){$(nRow).addClass('visited');};if(aData[18]=='L'){$(nRow).addClass('visitedLast')}})"),
         fnInfoCallback = I("(function( oSettings, iStart, iEnd, iMax, iTotal, sPre ) {Shiny.onInputChange('exploreTableiDSiDL',oSettings._iDisplayStart+'\t'+oSettings._iDisplayLength);return(sPre);})")
         #                                                          var p=0;
         #                                                          if(typeof Shiny.shinyapp.$inputValues.exploreTablePage != 'undefined'){p=Shiny.shinyapp.$inputValues.exploreTablePage};
         #                                                          var t=oSettings.oInstance;
         #                                                          if(Math.ceil(oSettings._iDisplayStart/oSettings._iDisplayLength)!=p){var tId=setTimeout(function(){t.fnPageChange(p,true)},100);}})")#,
         #            oTableTools = I("({'sSwfPath': 'TableTools-2.2.0/swf/copy_csv_xls_pdf.swf',
         #                              'aButtons': [{'sExtends':'text',
         #                                            'sButtonText':'Show all notes',
         #                                            'fnClick':function(button,config){
         #                                                           var oTable=getExploreTable();
         #                                                           $('#exploreTable tbody tr td img').each(function(){
         #                                                             var nTr=$(this).parents('tr')[0];
         #                                                             this.src='images/details_close.png';
         #                                                             oTable.fnOpen(nTr,fnFormatDetails(nTr),'details');
         #                                                           })
         #                                                      }
         #                                           },
         #                                           {'sExtends':'text',
         #                                            'sButtonText':'Hide all notes',
         #                                            'fnClick':function(button,config){
         #                                                           var oTable=getExploreTable();
         #                                                           $('#exploreTable tbody tr td img').each(function() {
         #                                                             var nTr=$(this).parents('tr')[0];
         #                                                             this.src='images/details_open.png';
         #                                                             oTable.fnClose(nTr);
         #                                                           })
         #                                                        }
         #                                            }
         #                                           ]
         #                               })")
    )
  }
  )
  
  observe({
    if (!is.null(input$excludeDeep) && (input$excludeDeep != isolate(values$lastExcludeDeep))) {
      values$lastExcludeDeep <- input$excludeDeep
      values$exploreTableiDS <- 0
    }
  })
  
  observe({
    if (!is.null(input$excludeExtra) && (input$excludeExtra != isolate(values$lastExcludeExtra))) {
      values$lastExcludeExtra <- input$excludeExtra
      values$exploreTableiDS <- 0
    }
  })
  
  observe({
    if (!is.null(input$excludeCommon) && (input$excludeCommon != isolate(values$lastExcludeCommon))) {
      values$lastExcludeCommon <- input$excludeCommon
      values$exploreTableiDS <- 0
    }
  })
  
  observe({
    if (!is.null(input$excludeIndels) && (input$excludeIndels != isolate(values$lastExcludeIndels))) {
      values$lastExcludeIndels <- input$excludeIndels
      values$exploreTableiDS <- 0
    }
  })
  
  
  ############################################
  # 3. Dynamic tabsetPanel for browsing info #
  #    contentType = "browseTabs"            #
  ############################################
  
  # Reactive function returning the info from CIGMA2 for the requested variation
  cigmainfo <- reactive({
    if (values$contentType=='exploreSingle' && !is.null(input$genetrans)){
      if (!nchar(input$genetrans)) {return()}
      gene<-unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))[1]
      transcript<-unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))[2]
      variation<-gsub(' [(].*$','',input$variationSel,perl=T)
    }else if (values$contentType=='browseTabs'){
      gene<-unlist(str_split(values$selectedVar,':'))[1]
      variation<-unlist(str_split(values$selectedVar,':'))[2]
    }else{
      return()
    }
    info(loggerServer,paste("in cigmainfo, gene=",gene," variation=",variation,sep=""))
    i<-NULL
    dataset<-NULL
    if (nchar(gene) && nchar(variation)) {
      mydb<-dbConnect(MySQL(),user='anonymous',dbname='cigma2',host=host)
      rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                 "select m.source,m.gene,m.transcript,rtranscript,hgvs_cdna,hgvs_prot,hgvs_prot_code1,
                m.altname,varLocation,r.name,varType,codingEffect,hg19_chr,hg19_pos,rsID,
                m.cdna_pos,offset,ntwt,ntmut,codon,aawt,aamut,c.firstorlast3,m.flags,pubmed,google_search_new,google_search_old
           from main_stable m 
      left join cdna2genomic c on m.gene=c.gene and m.cdna_pos=c.cdna_pos and m.varLocation='exon'
      left join capparegions r on m.gene=r.gene and m.hg19_chr=r.chr and m.hg19_pos between r.exonstart and r.exonend
          where m.gene='",gene,"' and m.hgvs_cdna='",variation,"'"))
      dataset<-fetch(rs,-1)
      if (!is.null(dataset) && dim(dataset)[1]>0) {
        i$rsid<-dataset[,'rsID']
        i$google_search_old<-dataset[,dim(dataset)[2]]
        dataset<-dataset[,-dim(dataset)[2]]
        i$google_search_new<-dataset[,dim(dataset)[2]]
        dataset<-dataset[,-dim(dataset)[2]]
        i$pubmed<-dataset[,dim(dataset)[2]]
        dataset<-dataset[,-dim(dataset)[2]]
        i$summary<-as.data.frame(t(dataset),stringsAsFactors=F)
        rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                   "SELECT type,IF(aa_from=aa_to,CONCAT(description,' [',aa_from,']'),CONCAT(description,' [',aa_from,'-',aa_to,']')) description 
           from main_stable m 
           join swissprot s on m.gene=s.gene and m.codon between s.aa_from and s.aa_to
          where m.gene='",gene,"' and m.hgvs_cdna='",variation,"'"))
        dataset<-fetch(rs,-1)
        if (!is.null(dataset) && dim(dataset)[1]>0) {
          i$swissprot<-as.data.frame(dataset,stringsAsFactors=F)
        }
        varsources <-strsplit(i$summary['source',],' ')[[1]]
        for (source in values$sources) {
          table <- tolower(source)
          if (source=='EVS') table <- 'esp_cappa'
          if (source=='ICR') table <- 'icr_controls'
          if (source=='HGMD') table <- 'hgmd_cappa'
          prefix <- table
          prefix <- sub('_cappa$','',prefix)
          prefix <- paste0(prefix,'_')
          if(source %in% varsources) {
            dataset<-NULL
            rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                       "select * from ",table," where gene='",gene,"' and hgvs_cdna='",variation,"'"))
            dataset<-fetch(rs,-1)
            dataset<-as.data.frame(t(dataset),stringsAsFactors=F)
            selectedrows<-grep(prefix,rownames(dataset))
            i[[tolower(source)]]<-as.data.frame(row.names=sub(paste0('^',prefix),'',rownames(dataset)[selectedrows]),
                                                dataset[selectedrows,],stringsAsFactors=F)
            if (source == 'TGP') {
              rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                         "select tgp_pop,tgp_AF,tgp_alleleCount,tgp_alleleTotal,tgp_genotypeCount,tgp_genotypePopSize from tgp_pop where tgp_rsID='",i$tgp["rsID",1],"'"))
              dataset<-fetch(rs,-1)
              i$tgp_pop<-as.data.frame(dataset,stringsAsFactors=F)
              colnames(i$tgp_pop)<-sub('^tgp_','',colnames(dataset))
              i$tgp_pop<-t(i$tgp_pop[order(i$tgp_pop$genotypePopSize,decreasing=T),])
            }
            if (source == 'TGP_PHASE3') {
              rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                         "select tgp_phase3_pop,tgp_phase3_AF,tgp_phase3_alleleCount,tgp_phase3_alleleTotal,tgp_phase3_genotypeCount,tgp_phase3_genotypePopSize from tgp_phase3_pop where tgp_phase3_id='",i$tgp_phase3["id",1],"'"))
              dataset<-fetch(rs,-1)
              i$tgp_phase3_pop<-as.data.frame(dataset,stringsAsFactors=F)
              colnames(i$tgp_phase3_pop)<-sub('^tgp_phase3_','',colnames(dataset))
              i$tgp_phase3_pop<-t(i$tgp_phase3_pop[order(i$tgp_phase3_pop$genotypePopSize,decreasing=T),])
            }
          }
        }
      }
      rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                 "select investigator_name,
             clinical_cigma_class,
             concat(date_format(creat_date,'%d/%m/%Y'),' (',
                    if(round(datediff(curdate(),creat_date)/7)<9,
                      concat(round(datediff(curdate(),creat_date)/7),' weeks ago'),
                      concat(round(datediff(curdate(),creat_date)/30.4166),' months ago')
                    ),')') date,
             notes
        from dogma_batchlog b 
       where b.gene='",gene,"' and b.variation='",variation,"'
       order by creat_date desc"));
      dataset<-fetch(rs,-1)
      if (!is.null(dataset) && dim(dataset)[1]>0) {
        i$notes<-as.data.frame(dataset,stringsAsFactors=F)
      }
      rs<-dbSendQuery(mydb,paste(sep="",collapse="",
                                 "select cigma_class, classification_justification
        from dogma_classification
       where gene='",gene,"' and hgvs_cdna='",variation,"' and version='",version,"'"));
      dataset<-fetch(rs,-1)
      if (!is.null(dataset) && dim(dataset)[1]>0) {
        i$baseline_class<-dataset[1,1]
        i$justification<-dataset[1,2]
      }
      dbDisconnect(mydb)
    }
    return(i)
  })
  
  # Output for the Summary panel
  output$summaryPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in summaryPanel")
    scholarlink1<-a(class="btn btn-primary",href=info$google_search_new,target="_blank","Google Scholar")
    scholarlink2<-a(class="btn btn-primary",href=info$google_search_old,target="_blank","Google Scholar (old nomen.)")
    search_new<-gsub('/scholar','/search',gsub('http://scholar.','http://www.',info$google_search_new))
    searchlink<-a(class="btn btn-primary",href=search_new,target="_blank","Google Search")
    pubmedlink <- ''
    if (!is.null(info) && !is.null(info$pubmed) && info$pubmed!=''){
      pubmedlink <- a(class="btn btn-primary",href=info$pubmed,target="_blank","HGMD Pubmed IDs")
    }
    dbsnplink <- ''
    clinvarlink <- ''
    if (!is.null(info) && !is.na(info$rsid) && info$rsid!=''){
      dbsnplink <- a(class="btn btn-primary",href=paste('http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',substr(info$rsid,3,nchar(info$rsid)),'#Diversity',sep=''),target="_blank","dbSNP")
      clinvarlink <- a(class="btn btn-primary",href=paste('http://www.ncbi.nlm.nih.gov/clinvar?term=',info$rsid,sep=''),target="_blank","ClinVar")
    }
    info.main<-as.list(info$summary[,1])
    names(info.main)<-rownames(info$summary)
    info.main[['name']] = gsub(info.main[['gene']],'',info.main[['name']])
    info.main[['name']] = gsub('on','on ',info.main[['name']])
    info.main[['name']] = gsub('of',' of ',info.main[['name']])
    firstorlast3<-fixempty(info.main[["firstorlast3"]])
    if (firstorlast3>0) firstorlast3<-paste0('+',as.character(firstorlast3))
    rmhclass<-''
    if (!is.null(info) && !is.null(info$rmh)) rmhclass<- info$rmh['class',1]
    sources<-gsub('ALAMUT', 'INSILICO', info.main[["source"]])
    sources<-gsub(' POLYPHEN2| MUTTASTER| CADD| SUSPECT','',sources)
    l<-list(tags$table(border="0",cellspacing="5",
                       tags$tr(tags$td(em(fields[["main"]][["gene"]]),             align="right"),tags$td(strong(info.main[["gene"]])),
                               tags$td(em(fields[["main"]][["hgvs_cdna"]]),        align="right"),tags$td(strong(info.main[["hgvs_cdna"]])),
                               tags$td(em(fields[["main"]][["varLocation"]]),      align="right"),tags$td(strong(info.main[["varLocation"]]))),
                       tags$tr(tags$td(em(fields[["main"]][["transcript"]]),       align="right"),tags$td(strong(info.main[["transcript"]])),
                               tags$td(em(fields[["main"]][["hgvs_prot"]]),        align="right"),tags$td(strong(fixempty(info.main[["hgvs_prot"]]))),
                               tags$td(em(fields[["main"]][["varType"]]),          align="right"),tags$td(strong(info.main[["varType"]]))),
                       tags$tr(tags$td(em(fields[["main"]][["rtranscript"]]),      align="right"),tags$td(strong(info.main[["rtranscript"]])),
                               tags$td(em(fields[["main"]][["hgvs_prot_code1"]]),  align="right"),tags$td(strong(fixempty(info.main[["hgvs_prot_code1"]]))),
                               tags$td(em(fields[["main"]][["codingEffect"]]),     align="right"),tags$td(strong(fixempty(info.main[["codingEffect"]])))),
                       tags$tr(tags$td(em(fields[["main"]][["hg19_chr"]]),         align="right"),tags$td(strong(info.main[["hg19_chr"]])),
                               tags$td(em(fields[["main"]][["altname"]]),          align="right"),tags$td(strong(fixempty(info.main[["altname"]]))),
                               tags$td(em(fields[["main"]][["cdna_pos"]]),         align="right"),tags$td(strong(info.main[["cdna_pos"]]))),
                       tags$tr(tags$td(em(fields[["main"]][["hg19_pos"]]),         align="right"),tags$td(strong(info.main[["hg19_pos"]])),
                               tags$td(em(fields[["main"]][["rsID"]]),             align="right"),tags$td(strong(fixempty(info.main[["rsID"]]))),
                               tags$td(em(fields[["main"]][["offset"]]),           align="right"),tags$td(strong(fixempty(info.main[["offset"]])))),
                       tags$tr(tags$td(em(fields[["main"]][["ntwt"]]),             align="right"),tags$td(strong(info.main[["ntwt"]])),
                               tags$td(em("Exon/intron:"),                         align="right"),tags$td(strong(info.main[["name"]])),
                               tags$td(em(fields[["main"]][["codon"]]),            align="right"),tags$td(strong(fixempty(info.main[["codon"]])))),
                       tags$tr(tags$td(em(fields[["main"]][["ntmut"]]),            align="right"),tags$td(strong(info.main[["ntmut"]])),
                               tags$td(),tags$td(),
                               tags$td(em(fields[["main"]][["firstorlast3"]]),     align="right"),tags$td(strong(firstorlast3)))
    ),
    br(),
    p(em('Data sources:'),strong(sources)),
    br(),
    p(em('Literature links:'),scholarlink1,scholarlink2,searchlink,pubmedlink),
    p(em('Database links:'),dbsnplink,clinvarlink),
    # p(em('Flags:'),strong(fixempty(info.main[["flags"]]))),
    # p(em('RMH class:'),strong(rmhclass)),
    br()
    )
    l1<-list()
    if (!is.null(info$hgmd) && info$hgmd['tag',]!='') {l1<-c(l1,list(p(em('HGMD classification:'),strong(info$hgmd['tag',]))))}
    if (!is.null(info$dmudb) && info$dmudb['interpretation',]!='') {l1<-c(l1,list(p(em('DMuDB classification:'),strong(info$dmudb['interpretation',]))))}
    if (!is.null(info$umd) && info$umd['significance',]!='') {l1<-c(l1,list(p(em('UMD classification:'),strong(info$umd['significance',]))))}
    l<-c(l,l1)
    return(as.character(tagList(l)))
  })
  
  # Output for the In silico predictions panel
  output$insilicoPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in insilicoPanel")
    l1a<-vector(mode="list")
    l2c<-vector(mode="list")
    l2s<-vector(mode="list")
    if (!is.null(info) && !is.null(info$alamut)) {
      fields.alamut<-names(fields[["alamut"]])  
      info.alamut<-as.list(info$alamut[,1])
      names(info.alamut)<-rownames(info$alamut)
      fields.missense<-fields.alamut[which(grepl('AGVGD|SIFT|MAPP',fields.alamut))]
      fields.conserv<-fields.alamut[which(grepl('Orthos|Conserved|BLOSUM',fields.alamut))]
      fields.splicing<-fields.alamut[which(grepl('SS|MaxEnt|NNS',fields.alamut))]
      l1a<-lapply(fields.missense,function(x){
        v<-info.alamut[[x]]
        if(v!='') {
          return(tags$tr(class=ifelse(x %in% c('AGVGDclass','SIFTprediction','MAPPprediction'),'highlight','normal'),tags$td(em(fields[["alamut"]][[x]]),align="right"),tags$td(strong(v))))
        }else{
          return('')
        }
      }
      )
      l2c<-lapply(fields.conserv,function(x){
        v<-info.alamut[[x]]
        if(v!='') {
          return(tags$tr(tags$td(em(fields[["alamut"]][[x]]),align="right"),tags$td(strong(v))))
        }else{
          return('')
        }
      }
      )
      l2s<-lapply(fields.splicing,function(x){
        v<-info.alamut[[x]]
        if(v!='') {
          return(tags$tr(tags$td(em(fields[["alamut"]][[x]]),align="right"),tags$td(strong(v))))
        }else{
          return('')
        }
      }
      )
    }
    l1m<-vector(mode="list")
    if (!is.null(info) && !is.null(info$muttaster)) {
      fields.muttaster<-names(fields[["muttaster"]])  
      info.muttaster<-as.list(info$muttaster[,1])
      names(info.muttaster)<-rownames(info$muttaster)
      l1m<-lapply(fields.muttaster,function(x){
        return(tags$tr(class=ifelse(x == 'prediction','highlight','normal'),tags$td(em(fields[["muttaster"]][[x]]),align="right"),tags$td(strong(info.muttaster[[x]]))))
      }
      )
    }
    l1p<-vector(mode="list")
    if (!is.null(info) && !is.null(info$polyphen2)) {
      fields.polyphen2<-names(fields[["polyphen2"]])  
      info.polyphen2<-as.list(info$polyphen2[,1])
      names(info.polyphen2)<-rownames(info$polyphen2)
      l1p<-lapply(fields.polyphen2,function(x){
        return(tags$tr(class=ifelse(x == 'hvar_prediction','highlight','normal'),tags$td(em(fields[["polyphen2"]][[x]]),align="right"),tags$td(strong(info.polyphen2[[x]]))))
      }
      )
    }
    l1c<-vector(mode="list")
    if (!is.null(info) && !is.null(info$cadd)) {
      fields.cadd<-names(fields[["cadd"]])  
      info.cadd<-as.list(info$cadd[,1])
      names(info.cadd)<-rownames(info$cadd)
      l1c<-lapply(fields.cadd,function(x){
        return(tags$tr(class=ifelse(x == 'Cscore','highlight','normal'),tags$td(em(fields[["cadd"]][[x]]),align="right"),tags$td(strong(info.cadd[[x]]))))
      }
      )
    }
    l1s<-vector(mode="list")
    if (!is.null(info) && !is.null(info$suspect)) {
      fields.suspect<-names(fields[["suspect"]])  
      info.suspect<-as.list(info$suspect[,1])
      names(info.suspect)<-rownames(info$suspect)
      l1s<-lapply(fields.suspect,function(x){
        return(tags$tr(class=ifelse(x == 'score','highlight','normal'),tags$td(em(fields[["suspect"]][[x]]),align="right"),tags$td(strong(info.suspect[[x]]))))
      }
      )
    }
    l1<-list(h4('Missense/nonsense predictions'),
             tags$table(border="0",cellspacing="5",l1a,l1m,l1p,l1c,l1s)
    )
    
    l2<-vector(mode="list")
    if (!is.null(info) && !is.null(info$swissprot)) {
      l2<-list(h4('Swissprot features'),
               tags$table(border="0",cellspacing="5",
                          tags$tr(tags$td(em('Type')),tags$td(em('Description'))),
                          lapply(1:dim(info$swissprot)[1],function(x){
                            return(tags$tr(tags$td(info$swissprot[x,1]),tags$td(info$swissprot[x,2])))
                          }
                          )
               )
      )
    }
    l2<-c(l2,list(h4('Conservation'),tags$table(border="0",cellspacing="5",l2c),
                  h4('Splicing predictions'),tags$table(border="0",cellspacing="5",l2s)
    )
    )
    return(as.character(tagList(fluidRow(column(width=8,l1),column(width=4,l2)))))
  })
  
  # Output for the Frequency panel
  output$frequencyPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in frequencyPanel")
    l<-vector(mode="list")
    if (!is.null(info) && !is.null(info$evs)) {
      fields.evs<-names(fields[["esp"]])
      info.evs<-as.list(info$evs[,1])
      names(info.evs)<-rownames(info$evs)
      info.main<-as.list(info$summary[,1])
      names(info.main)<-rownames(info$summary)
      l<-c(l,list(h4(a('Exome Sequencing Project (ESP) / Exome Variant Server (EVS)',
                       href=paste0('http://evs.gs.washington.edu/EVS/PopStatsServlet?searchBy=chromosome&chromosome=',
                                   info.main[["hg19_chr"]],'&chromoStart=',info.main[["hg19_pos"]],
                                   '&chromoEnd=',info.main[["hg19_pos"]],'&x=0&y=0'),target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.evs,function(x){
                               v<-info.evs[[x]]
                               if (grepl('MAF',x)) v<-freqprint(v)
                               return(tags$tr(tags$td(em(fields[["esp"]][[x]]),align="right"),tags$td(strong(v))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$tgp)) {
      fields.tgp<-names(fields[["tgp"]])  
      info.tgp<-as.list(info$tgp[,1])
      names(info.tgp)<-rownames(info$tgp)
      fields.tgp_pop<-names(fields[["tgp_pop"]])
      l<-c(l,list(h4(a('1000 Genome project (TGP / 1KG)',
                       href=paste0('http://browser.1000genomes.org/Homo_sapiens/Variation/Population?db=core;v=',info$rsid,';vdb=variation'),
                       target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.tgp,function(x){
                               v<-info.tgp[[x]]
                               if (grepl('AF',x)) v<-freqprint(v)
                               return(tags$tr(tags$td(em(fields[["tgp"]][[x]]),align="right"),tags$td(strong(v))))
                             }
                             )
                  ),
                  tags$table(border="1",
                             lapply(fields.tgp_pop,function(x){
                               return(tags$tr(tags$td(em(fields[["tgp_pop"]][[x]]),align="right"),
                                              lapply(info$tgp_pop[x,],function(y,fld=x){
                                                v<-y
                                                if (grepl('AF',fld)) v<-freqprint(v)
                                                if (fld=='pop'){
                                                  return(tags$td(strong(v),align="right",title=tp[[v]]))
                                                }else{
                                                  return(tags$td(strong(v),align="right"))
                                                }
                                              }
                                              )
                               )
                               )
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$tgp_phase3)) {
      fields.tgp_phase3<-names(fields[["tgp_phase3"]])  
      info.tgp_phase3<-as.list(info$tgp_phase3[,1])
      names(info.tgp_phase3)<-rownames(info$tgp_phase3)
      fields.tgp_phase3_pop<-names(fields[["tgp_phase3_pop"]])
      l<-c(l,list(h4('1000 Genome project - phase3 (TGP-phase3 / 2.5KG)'),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.tgp_phase3,function(x){
                               v<-info.tgp_phase3[[x]]
                               if (grepl('AF',x)) v<-freqprint(v)
                               return(tags$tr(tags$td(em(fields[["tgp_phase3"]][[x]]),align="right"),tags$td(strong(v))))
                             }
                             )
                  ),
                  tags$table(border="1",
                             lapply(fields.tgp_phase3_pop,function(x){
                               return(tags$tr(tags$td(em(fields[["tgp_phase3_pop"]][[x]]),align="right"),
                                              lapply(info$tgp_phase3_pop[x,],function(y,fld=x){
                                                v<-y
                                                if (grepl('AF',fld)) v<-freqprint(v)
                                                if (fld=='pop'){
                                                  return(tags$td(strong(v),align="right",title=tp[[v]]))
                                                }else{
                                                  return(tags$td(strong(v),align="right"))
                                                }
                                              }
                                              )
                               )
                               )
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$icr)) {
      fields.icr<-names(fields[["icr"]])  
      info.icr<-as.list(info$icr[,1])
      names(info.icr)<-rownames(info$icr)
      l<-c(l,list(h4('1958 Birth Cohort (ICR)'),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.icr,function(x){
                               return(tags$tr(tags$td(em(fields[["icr"]][[x]]),align="right"),tags$td(strong(info.icr[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$bic)) {
      fields.bic<-names(fields[["bic"]])  
      info.bic<-as.list(info$bic[,1])
      names(info.bic)<-rownames(info$bic)
      l<-c(l,list(h4('BIC'),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.bic,function(x){
                               return(tags$tr(tags$td(em(fields[["bic"]][[x]]),align="right"),tags$td(strong(info.bic[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$lovd)) {
      fields.lovd<-names(fields[["lovd"]])
      fields.lovd<-fields.lovd[-which(fields.lovd=='insight_class')]
      info.lovd<-as.list(info$lovd[,1])
      names(info.lovd)<-rownames(info$lovd)
      gene<-info$summary['gene',1]
      variant<-info$summary['hgvs_cdna',1]
      lovd_link<-'LOVD'
      if (grepl('BRCA',gene)){
        lovd_link <- a(href=paste0('http://chromium.liacs.nl/LOVD2/cancer/variants.php?select_db=',gene,'&action=search_unique&search_Variant%2FDNA=',variant),'LOVD',target="_blank") 
      }
      if (gene %in% c('MLH1','MSH2','MSH6','PMS2')){
        lovd_link <- a(href=paste0('http://chromium.liacs.nl/LOVD2/colon_cancer/variants.php?select_db=',gene,'&action=search_unique&search_Variant%2FDNA=',variant),'LOVD',target="_blank")
      }
      l<-c(l,list(h4(lovd_link),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.lovd,function(x){
                               return(tags$tr(tags$td(em(fields[["lovd"]][[x]]),align="right"),tags$td(strong(info.lovd[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$dmudb)) {
      fields.dmudb<-names(fields[["dmudb"]])  
      info.dmudb<-as.list(info$dmudb[,1])
      names(info.dmudb)<-rownames(info$dmudb)
      l<-c(l,list(h4('DMuDB'),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.dmudb,function(x){
                               return(tags$tr(tags$td(em(fields[["dmudb"]][[x]]),align="right"),tags$td(strong(info.dmudb[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$umd)) {
      fields.umd<-names(fields[["umd"]])  
      info.umd<-as.list(info$umd[,1])
      names(info.umd)<-rownames(info$umd)
      l<-c(l,list(h4('UMD'),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.umd,function(x){
                               return(tags$tr(tags$td(em(fields[["umd"]][[x]]),align="right"),tags$td(strong(info.umd[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    return(as.character(tagList(l)))
  })
  
  # Output for the Genetic Epidemiology panel
  output$geneticepiPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in geneticepiPanel")
    l<-vector(mode="list")
    if (!is.null(info) && !is.null(info$easton)) {
      fields.easton<-names(fields[["easton"]])
      info.easton<-as.list(info$easton[,1])
      names(info.easton)<-rownames(info$easton)
      l<-c(l,list(h4(a('Easton et al, Am J Hum Genet. 2007',href='http://www.ncbi.nlm.nih.gov/pubmed/17924331',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.easton,function(x){
                               return(tags$tr(tags$td(em(fields[["easton"]][[x]]),align="right"),tags$td(strong(info.easton[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$lindor)) {
      fields.lindor<-names(fields[["lindor"]])
      info.lindor<-as.list(info$lindor[,1])
      names(info.lindor)<-rownames(info$lindor)
      l<-c(l,list(h4(a('Lindor et al, Hum Mutat. 2011',href='http://www.ncbi.nlm.nih.gov/pubmed/21990134',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.lindor,function(x){
                               v<-info.lindor[[x]]
                               if (grepl('reference',x) && info$lindor['pubmed',1]!='') v<-a(v,href=paste0('http://www.ncbi.nlm.nih.gov/pubmed/',info$lindor['pubmed',1]),target='_blank')
                               return(tags$tr(tags$td(em(fields[["lindor"]][[x]]),align="right"),tags$td(strong(v))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$iarc)) {
      fields.iarc<-names(fields[["iarc"]])
      info.iarc<-as.list(info$iarc[,1])
      names(info.iarc)<-rownames(info$iarc)
      l<-c(l,list(h4(a('IARC: Vallee et al, Hum Mutat. 2012',href='http://www.ncbi.nlm.nih.gov/pubmed/21990165',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.iarc,function(x){
                               return(tags$tr(tags$td(em(fields[["iarc"]][[x]]),align="right"),tags$td(strong(info.iarc[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$lovd) && info$lovd['insight_class',]!='') {
      gene<-info$summary['gene',1]
      variant<-info$summary['hgvs_cdna',1]
      l<-c(l,list(h4(a('InSiGHT: Thompson et al, Nat Genet. 2014',href='http://www.ncbi.nlm.nih.gov/pubmed/24362816',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             tags$tr(tags$td(em('InSiGHT class:'),align="right"),tags$td(a(strong(info$lovd['insight_class',]),href=paste0('https://googledrive.com/host/0B8HVsr5izQxJUi1XTzEtWFlRc00/index.html?gene=',gene,'&protein=&variant=',variant),target='_blank')))
                  )
      )
      )
    }
    return(as.character(tagList(l)))
  })
  
  # Output for the Splicing Analysis panel
  output$splicingPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in splicingPanel")
    l<-vector(mode="list")
    if (!is.null(info) && !is.null(info$houdayer)) {
      fields.houdayer<-names(fields[["houdayer"]])
      info.houdayer<-as.list(info$houdayer[,1])
      names(info.houdayer)<-rownames(info$houdayer)
      l<-c(l,list(h4(a('Houdayer et al, Hum Mutat. 2012',href='http://www.ncbi.nlm.nih.gov/pubmed/22505045',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.houdayer,function(x){
                               return(tags$tr(tags$td(em(fields[["houdayer"]][[x]]),align="right"),tags$td(strong(info.houdayer[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$walker)) {
      fields.walker<-names(fields[["walker"]])
      l<-c(l,list(h4(a('Walker et al, Hum Mutat. 2013',href='http://www.ncbi.nlm.nih.gov/pubmed/23893897',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.walker,function(x){
                               return(tags$tr(tags$td(em(fields[["walker"]][[x]]),align="right"),
                                              lapply(info$walker[x,],function(y,fld=x){
                                                v<-y
                                                if (grepl('pubmed',fld)) v<-a(v,href=paste0('http://www.ncbi.nlm.nih.gov/pubmed/',v),target='_blank')
                                                return(tags$td(strong(v),align="right"))
                                              }
                                              )
                               )
                               )
                             }
                             )
                  )
      )
      )
    }
    return(as.character(tagList(l)))
  })
  
  # Output for the Functional Analysis panel
  output$functionalPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in functionalPanel")
    l<-vector(mode="list")
    if (!is.null(info) && !is.null(info$guidugli)) {
      fields.guidugli<-names(fields[["guidugli"]])
      info.guidugli<-as.list(info$guidugli[,1])
      names(info.guidugli)<-rownames(info$guidugli)
      l<-c(l,list(h4(a('Guidugli et al, Cancer Res. 2013',href='http://www.ncbi.nlm.nih.gov/pubmed/23108138',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.guidugli,function(x){
                               return(tags$tr(tags$td(em(fields[["guidugli"]][[x]]),align="right"),tags$td(strong(info.guidugli[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$drost)) {
      fields.drost<-names(fields[["drost"]])
      info.drost<-as.list(info$drost[,1])
      names(info.drost)<-rownames(info$drost)
      if (grepl('MSH',info$summary['gene',1])) {
        title<-a('Drost et al, Hum Mutat. 2012',href='http://www.ncbi.nlm.nih.gov/pubmed/22102614',target='_blank')
      }else{
        title<-a('Drost et al, Hum Mutat. 2010',href='http://www.ncbi.nlm.nih.gov/pubmed/20020535',target='_blank')  
      }
      l<-c(l,list(h4(title),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.drost,function(x){
                               return(tags$tr(tags$td(em(fields[["drost"]][[x]]),align="right"),
                                              lapply(info.drost[[x]],function(y,fld=x){
                                                v<-y
                                                if (grepl('reference|ficient',fld)) v<-gsub('([^(,]+)\\(PMID:(\\d+)\\)','<A HREF="http://www.ncbi.nlm.nih.gov/pubmed/\\2" target="_blank">\\1</A>',y,perl=T)
                                                return(tags$td(strong(HTML(v))))
                                              }
                                              )
                               )
                               )
                             }
                             )
                  )
      )
      )
    }
    if (!is.null(info) && !is.null(info$bouwman)) {
      fields.bouwman<-names(fields[["bouwman"]])
      info.bouwman<-as.list(info$bouwman[,1])
      names(info.bouwman)<-rownames(info$bouwman)
      l<-c(l,list(h4(a('Bouwman et al, Cancer Discov. 2013',href='http://www.ncbi.nlm.nih.gov/pubmed/23867111',target='_blank')),
                  tags$table(border="0",cellspacing="5",
                             lapply(fields.bouwman,function(x){
                               return(tags$tr(tags$td(em(fields[["bouwman"]][[x]]),align="right"),tags$td(strong(info.bouwman[[x]]))))
                             }
                             )
                  )
      )
      )
    }
    
    return(as.character(tagList(l)))
  })
  
  output$notesPanel <- renderText ({
    info<-cigmainfo()
    info(loggerServer,"in notesPanel")
    l<-list();
    notes_new <- NULL
    if (values$contentType=='browseTabs') {
      notes_new <- values$exploreTable[which(grepl(values$selectedVar,values$exploreTable$variation,fixed=TRUE)),'notes_new']
      if (notes_new[1] != '') {
        nn <- unlist(str_split(notes_new[1],"\t"))
        ntext <- gsub('PMID[:]?(\\s+)?(\\d+)','<A HREF="http://www.ncbi.nlm.nih.gov/pubmed?term=\\2" target="_blank">PMID:\\2</A>',nn[4],perl=T)
        l<-c(l,
             list(div(class='highlight',
                      fluidRow(column(width=4,em('Curated class:')),column(width=8,strong(nn[2]))),
                      fluidRow(column(width=4,em('Date: ')),column(width=8,nn[3])),
                      fluidRow(column(width=4,em('Name: ')),column(width=8,nn[1]))),
                  fluidRow(div(class="details",style="margin-bottom:20px",lapply(unlist(str_split(ntext,'<br>')),function(x){div(HTML(x))})))
             )
        )
      }
    }
    if (!is.null(info) && !is.null(info$notes) && !is.na(info$notes)) {
      l<-c(l,apply(info$notes,1,function(x){
        ntext <- gsub('PMID:(\\s+)?(\\d+)','<A HREF="http://www.ncbi.nlm.nih.gov/pubmed?term=\\2" target="_blank">PMID:\\2</A>',x[4],perl=T)
        list(div(class='highlight',
                 fluidRow(column(width=4,em('Curated class: ')),column(width=8,strong(x[2]))),
                 fluidRow(column(width=4,em('Date: ')),column(width=8,x[3])),
                 fluidRow(column(width=4,em('Name: ')),column(width=8,x[1]))),
             fluidRow(div(class="details",style="margin-bottom:20px",em('Notes:'),lapply(unlist(str_split(ntext,'<br>')),function(x){div(HTML(x))})))
        )
      }
      )
      )
    }
    l<-c(l,list(div(class='highlight',
                    fluidRow(column(width=4,em('Automated class:')),column(width=8,strong(info$baseline_class))),
                    fluidRow(column(width=4,em('Justification:')),column(width=8,strong(info$justification)))
    )
    )
    )
    return(as.character(tagList(l)))
  })  
  
  ###################
  # Downloaded data #
  ###################
  output$downloadVariantReport = downloadHandler(
    filename =  function() {info<-cigmainfo()
                            paste0(info$summary['gene',1],'_',info$summary['hgvs_cdna',1],'_variantReport.pdf') },
    content = function(file) {
      debug(loggerServer, 'VariantReport: getting cigmainfo')
      info<-cigmainfo()
      debug(loggerServer, '               getting VariantReport')
      variantReport<-getVariantReport()
      debug(loggerServer, '               getting VariantNotes')
      variantNotes<-getVariantNotes()
      debug(loggerServer, '          done getting VariantNotes')
      savewd<-getwd()
      debug(loggerServer, paste('current directory=',savewd))
      debug(loggerServer, '          done getting VariantNotes')
      debug(loggerServer, 'dumping variables')
      dump(c('info','variantReport','variantNotes'),file=paste(savewd,'reports','dump.txt',sep='/'))
      debug(loggerServer, 'dumping done')
#      setwd('reports') # knitr can only write in folder reports
#      tryCatch({
        debug(loggerServer, '               creating .tex from .Rnw')
        knit('reports/cavada_variant_report.Rnw')
        debug(loggerServer, '               creating .pdf from .tex')
        system('/usr/bin/pdflatex cavada_variant_report.tex')
  
        debug(loggerServer, paste('               renaming cavada_variant_report.pdf to',file))
        file.rename('cavada_variant_report.pdf', file) # move pdf to file for downloading
#      },finally = {
#        setwd(savewd)
        debug(loggerServer, 'done creating .pdf from .Rnw')
#      })
    },
    
    contentType = 'application/pdf'
  )
  
  output$downloadAnnotatedBatch <- downloadHandler(
    filename = function() { paste0(values$iname,'_',values$sessionid,'.curated.csv') },
    content = function(file) {write.csv(getAnnotatedBatch(), file, row.names=F, quote=T, eol="\r\n")}
  )
  
  getAnnotatedBatch <- reactive({
    if (values$inputFileType=='simple') {
      v<-data.frame(sampleID=values$exploreTable[,'sample_id'],
                    gene=sapply(values$exploreTable[,'variation'],function(x){return(unlist(str_split(x,':'))[1])}),
                    hgvs_cdna=sapply(values$exploreTable[,'variation'],function(x){return(unlist(str_split(unlist(str_split(x,':'))[2],' '))[1])}),
                    hgvs_protein=sapply(values$exploreTable[,'variation'],function(x){return(gsub('[()]','',unlist(str_split(unlist(str_split(x,':'))[2],' '))[2]))}),
                    varType=values$exploreTable[,'varType'],
                    varLocation=values$exploreTable[,'varLocation'],
                    codingEffect=values$exploreTable[,'codingEffect'],
                    automated_class=values$exploreTable[,'working_cigma_class'],
                    classification_justification=values$exploreTable[,'classification_justification'],
                    curated_class=values$exploreTable[,'clinical_cigma_class'])
      s<-lapply(1:dim(v)[1],
                function(y){
                  return(sapply(names(fan),
                                function(x,yy=y){
                                  table<-x
                                  if(table=='esp') table<-'esp_cappa'
                                  if(table=='hgmd') table<-'hgmd_cappa'
                                  prefix<-x
                                  if(table=='tgp_pop') {
                                    prefix<-'tgp'
                                    table<-'tgp_pop_main'
                                  }
                                  if(table=='tgp_phase3_pop') {
                                    prefix<-'tgp_phase3'
                                    table<-'tgp_phase3_pop_main'
                                  }
                                  fields <- paste(paste0(prefix,'_',fan[[x]]),collapse=',')
                                  if(table=='walker') fields <- paste0('GROUP_CONCAT(',fields,') ',fields)
                                  select_statement<-paste0("select ",fields," from ",table," where gene='",v[yy,2],"' and hgvs_cdna='",v[yy,3],"'")
                                  return(select_statement)
                                }
                  )
                  )
                }
      )
      mydb<-dbConnect(MySQL(),user='anonymous',dbname='cigma2',host=host)
      extra<-lapply(1:length(s),
                    function(y) {
                      return(
                        unlist(sapply(1:length(s[[y]]),function(x){rs<-dbSendQuery(mydb,s[[y]][x]);dataset<-fetch(rs,-1);return(as.vector(dataset))}))
                      )
                    }
      )
      dbDisconnect(mydb)
      for(i in 1:length(extra)) {
        if(length(extra[[i]])>0) {
          for (j in 1:length(extra[[i]])) {
            v[i,names(extra[[i]])[j]]<-extra[[i]][j]
          }  
        }
      }
      for (j in 1:dim(v)[2]) {
        if (length(which(is.na(v[,j])))>0) {
          v[which(is.na(v[,j])),j]<-''
        }
      }
      colnames(v)=c('Sample ID',
                    'Gene',
                    'HGVS cDNA',
                    'HGVS protein',
                    'Variant type',
                    'Variant location',
                    'Variant coding effect',
                    'Automated class',
                    'Class justification',
                    'Curated class',
                    fa[colnames(v)[11:length(colnames(v))],3])
      ordered_fields<-fa[order(fa$analysis_output_order),'field_label']
      return(v[,c(colnames(v)[1:10],ordered_fields[which(ordered_fields %in% colnames(v))])])
    }
    if (values$inputFileType=='clinical') {
      v<-data.frame(Worksheet=values$exploreTable[,'worksheet'],
                    InvID=values$exploreTable[,'inv_id'],
                    DNA.No=values$exploreTable[,'sample_id'],
                    Gene=sapply(values$exploreTable[,'variation'],function(x){return(unlist(str_split(x,':'))[1])}),
                    Nomenclature=sapply(values$exploreTable[,'variation'],function(x){return(unlist(str_split(x,':'))[2])}),
                    Het.Hom=values$exploreTable[,'zygosity'],
                    Class=values$exploreTable[,'clinical_cigma_class'],
                    Exon=values$exploreTable[,'exon'])
      return(v)
    }
  })
  
  getVariantReport <- reactive({
    if (values$contentType=='exploreSingle' && !is.null(input$genetrans)){
      if (!nchar(input$genetrans)) {return()}
      gene<-unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))[1]
      transcript<-unlist(str_extract_all(input$genetrans,pattern="[A-Z0-9]+"))[2]
      variation<-gsub(' [(].*$','',input$variationSel,perl=T)
    }else if (values$contentType=='browseTabs'){
      gene<-unlist(str_split(values$selectedVar,':'))[1]
      variation<-unlist(str_split(values$selectedVar,':'))[2]
    }else{
      return()
    }
    debug(loggerServer,paste("in getVariantReport, gene=",gene," variation=",variation,sep=""))
    ordered_fields<-rownames(fv)[order(fv$variant_report_order)]
    info<-cigmainfo()
    debug(loggerServer,'   got cigmainfo')
    info$summary['name',1] = gsub(info$summary['gene',1],'',info$summary['name',1])
    info$summary['name',1] = gsub('on','on ',info$summary['name',1])
    info$summary['name',1] = gsub('of',' of ',info$summary['name',1])
    variantReport<-data.frame(
      Description=c(
        'Gene',
        'Ensembl transcript',
        'RefSeq transcript',
        'Chromosome',
        'Genomic position',
        'Reference allele',
        'Variant allele',
        'HGVS cDNA',
        'HGVS protein',
        'HGVS protein (1 letter)',
        'Alternative notation',
        'rsID',
        'Variant location',
        'Exon/intron number',
        'Variant type',
        'Variant coding effect',
        'cDNA position',
        'Codon number',
        'Intronic position',
        'First/last 3 bases of exon'
      ),
      Value=c(
        gene,
        info$summary['transcript',1],
        info$summary['rtranscript',1],
        info$summary['hg19_chr',1],
        info$summary['hg19_pos',1],
        info$summary['ntwt',1],
        info$summary['ntmut',1],
        variation,
        info$summary['hgvs_prot',1],
        info$summary['hgvs_prot_code1',1],
        info$summary['altname',1],
        info$summary['rsID',1],
        info$summary['varLocation',1],
        info$summary['name',1],
        info$summary['varType',1],
        info$summary['codingEffect',1],
        info$summary['cdna_pos',1],
        info$summary['codon',1],
        info$summary['offset',1],
        info$summary['firstorlast3',1]
      ),
      row.names=c(
        'main_gene',
        'main_transcript',
        'main_rtranscript',
        'main_hg19_chr',
        'main_hg19_pos',
        'main_ntwt',
        'main_ntmut',
        'main_hgvs_cdna',
        'main_hgvs_prot',
        'main_hgvs_prot_code1',
        'main_altname',
        'main_rsID',
        'main_varLocation',
        'main_name',
        'main_varType',
        'main_codingEffect',
        'main_cdna_pos',
        'main_codon',
        'main_offset',
        'main_firstorlast3'
      )
    )
    # hack to insert the exon/intron number after varLocation
    ordered_fields<-c(ordered_fields[1:13],'main_name',ordered_fields[14:length(ordered_fields)])
    if (!is.null(info) && !is.null(info$swissprot)) {
      swissprots<-lapply(1:dim(info$swissprot)[1],function(x){
        return(c(paste0('Swissprot feature',x,' - ',info$swissprot[x,1]),info$swissprot[x,2]))
      })
      swissprot_fields<-NULL
      for(i in 1:length(swissprots)) {
        variantReport<-rbind(variantReport,data.frame(Description=swissprots[[i]][1],Value=swissprots[[i]][2],row.names=paste0('swissprot_feature_',i)))
        swissprot_fields<-c(swissprot_fields,paste0('swissprot_feature_',i))
      }
      # hack to insert the swissprot features after firstorlast3
      ordered_fields<-c(ordered_fields[1:17],swissprot_fields,ordered_fields[18:length(ordered_fields)])
    }
    
    t<-names(fvn)
    t<-t[-which(t=='main')]
    s<-lapply(t,
              function(x){
                table<-x
                if(table=='esp') table<-'esp_cappa'
                if(table=='hgmd') table<-'hgmd_cappa'
                prefix<-x
                if(table=='icr') {
                  table<-'icr_controls'
                  prefix<-'icr_controls'
                }  
                if(table=='tgp_pop') {
                  prefix<-'tgp'
                  table<-'tgp_pop_main'
                }
                if(table=='tgp_phase3_pop') {
                  prefix<-'tgp_phase3'
                  table<-'tgp_phase3_pop_main'
                }
                fields <- paste(paste0(prefix,'_',fvn[[x]]),collapse=',')
                if(table=='walker') fields <- paste0('GROUP_CONCAT(',fields,') ',fields)
                select_statement<-paste0("select ",fields," from ",table," where gene='",gene,"' and hgvs_cdna='",variation,"'")
                return(select_statement)
              }
    )
    s<-setNames(s,names(t))
    mydb<-dbConnect(MySQL(),user='anonymous',dbname='cigma2',host=host)
    extra<-lapply(1:length(s),function(x){
                                          debug(loggerServer,paste('trying to run SQL statement:',s[[x]],sep='\n'))
                                          rs<-dbSendQuery(mydb,s[[x]]);
                                          dataset<-fetch(rs,-1);
                                          return(as.vector(dataset))
                                         }
                 )
    extra<-setNames(extra,names(t))
    dbDisconnect(mydb)
    debug(loggerServer,'adding extra fields')
    for(i in 1:length(extra)) {
      if(length(extra[[i]])>0) {
        for (j in 1:length(extra[[i]])) {
          field_label <- fv[names(extra[[i]])[j],'field_label']
          if (!is.na(field_label) && !is.na(extra[[i]][1,j])) {
            variantReport<-rbind(variantReport,data.frame(Description=field_label,Value=extra[[i]][1,j],row.names=names(extra[[i]])[j]))
          }
        }
      }
    }
    debug(loggerServer,'done adding extra fields')
    if (!(gene %in% c('MSH2','MSH6','MLH1','PMS2')) && ('lovd_insight_class' %in% row.names(variantReport))){
      variantReport<-variantReport[-which(row.names(variantReport)=='lovd_insight_class'),]
    }
    variantReport<-variantReport[ordered_fields[which(ordered_fields %in% row.names(variantReport))],]
    variantReport<-variantReport[-which(is.na(variantReport$Value)),]
    variantReport[,1]<-paste0(variantReport[,1],':')
    variantReport1<-data.frame(Values=variantReport[,2],row.names=variantReport[,1])
    debug(loggerServer,'done fixing the variantReport')
    return(variantReport1)
  })
  
  sanitize <- function(str) {
    Sys.setlocale('LC_ALL','C')
    result <- str
    result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
    result <- gsub("$", "\\$", result, fixed = TRUE)
    result <- gsub(">", "$>$", result, fixed = TRUE)
    result <- gsub("<", "$<$", result, fixed = TRUE)
    result <- gsub("|", "$|$", result, fixed = TRUE)
    result <- gsub("{", "\\{", result, fixed = TRUE)
    result <- gsub("}", "\\}", result, fixed = TRUE)
    result <- gsub("%", "\\%", result, fixed = TRUE)
    result <- gsub("&", "\\&", result, fixed = TRUE)
    result <- gsub("_", "\\_", result, fixed = TRUE)
    result <- gsub("#", "\\#", result, fixed = TRUE)
    result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
    result <- gsub("~", "\\~{}", result, fixed = TRUE)
    result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$", result, fixed = TRUE)
    return(result)
  }
  
  getVariantNotes<-reactive({
    info<-cigmainfo()
    info(loggerServer,'in getVariantNotes')
    variantNotes<-data.frame()
    if (!is.null(info) && !is.null(info$notes) && !is.na(info$notes)) {
      mynotes<-info$notes[,c('date','notes')]
      debug(loggerServer,'sanitizing mynotes')
      mynotes$notes<-sanitize(mynotes$notes)
      debug(loggerServer,'adding PMID URLs')
      mynotes$notes<-gsub('PMID[:]?(\\s+)?(\\d+)','\\\\href{http://www.ncbi.nlm.nih.gov/pubmed?term=\\2}{PMID:\\2}',gsub('<br>','\\baselineskip',mynotes$notes))
      colnames(mynotes)<-c('Date','Notes')
      variantNotes<-rbind(variantNotes,mynotes)
    }
    debug(loggerServer,'done with notes')
    return(variantNotes)
  })
  #################
  # Other outputs #
  #################
  
  output$exploreTableVariation <- renderText({
    input$exploreTableVariation
  })
  
  # Output variation name (for the browseTab screen)
  output$selectedVar <- renderText({
    values$selectedVar
  })
  
  output$selectedVarBaseline <- renderText({
    ifelse(values$selectedVarBaseline!='',values$selectedVarBaseline,' ')
  })
  
  output$selectedSet <- renderText({
    input$batchid
  })
  
})
