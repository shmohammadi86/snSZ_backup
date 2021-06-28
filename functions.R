synapse.project.id = "syn22755334"
synapse.results.folder = "syn25910026"

Syn.datasets <- Folder(name = "datasets", parentId = synapse.results.folder)
Syn.tables   <- Folder(name = "tables", parentId = synapse.results.folder)
Syn.figures  <- Folder(name = "figures", parentId = synapse.results.folder)
Syn.input    <- Folder(name = "input", parentId = synapse.project.id)

# Functions to Pull/Push objects from/to Synapse
storeDataset = function( obj, name, description = NULL, dataset.path = "results/datasets"){
  file_name = paste0(name, ".rds")
  file_path = file.path(dataset.path, file_name)
  readr::write_rds( obj, file=file_path)
  
  Syn.datasets <- synStore(Syn.datasets)
  folder = Syn.datasets$properties$id
  if(is.null(description)) {
    description = file_name
  }	
  OBJ = File(file_path, name = description, parentId = folder)
  
  suppressWarnings( {suppressMessages( {out = synStore(OBJ)} )} )	
}

# Store Table
storeTable = function( DFs, name, description = NULL, tables.path = "results/tables"){
  if(is.null(names(DFs))) {
    names(DFs) = paste("Table", 1:length(DFs), sep = "")
  }
  wb <- openxlsx::createWorkbook()
  for(i in 1:length(DFs)) {
    n = names(DFs)[[i]] 
    openxlsx::addWorksheet(wb=wb, sheetName = n)
    openxlsx::writeData(wb, sheet = n, DFs[[i]]) 
  }
  file_name = paste0(name, ".xlsx")
  file_path = file.path(tables.path, file_name)
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  
  Syn.tables <- synStore(Syn.tables)
  folder = Syn.tables$properties$id
  if(is.null(description)) {
    description = file_name
  }
  OBJ = File(path = file_path, name = description, parentId = folder)
  
  suppressWarnings( {suppressMessages( {out = synStore(OBJ)} )} )	
}

# Store figure
storeFigure = function( gg, name, extension, description = NULL, width = 10, height = 8, figures.path = "results/figures"){
  if(extension == "png") {
    file_name = paste0(name, ".png")
    file_path = file.path(figures.path, file_name)
    png(file_path, width = width, height = height, units = "in", res = 300)
  } else if(extension == "pdf") {
    file_name = paste0(name, ".pdf")
    file_path = file.path(figures.path, file_name)
    pdf(file_path, width = width, height = height)
  }
  print(gg)
  dev.off()
  
  Syn.figures <- synStore(Syn.figures)
  folder = Syn.figures$properties$id
  if(is.null(description)) {
    description = file_name
  }  
  OBJ = File(path = file_path, name = description, parentId = folder)
  
  suppressWarnings( {suppressMessages( {out = synStore(OBJ)} )} )
}

# Load a dataset, either locally or fetch from Synapse
loadDataset = function( name, dataset.path = "results/datasets" ){
  file_name = paste0(name, ".rds")
  file_path = file.path(dataset.path, file_name)
  if(!file.exists(file_path)) {
    Syn.datasets <- synStore(Syn.datasets)
    folder = Syn.datasets$properties$id
    
    syn.files = synChildren(folder)
    mask = names(syn.files) == file_name
    if(sum(mask) == 0) {
      warning(sprintf("Did not find file %s", file_name))
      return()
    } else {
      syn.id = as.character(syn.files[which(mask == T)[[1]]])
      entity <- synGet(syn.id, downloadLocation = dataset.path)  	
    }
  }
  obj = readr::read_rds( file_path )
  
  return(obj)
}

# Load an "input" dataset, either locally or fetch from Synapse
loadInputDataset = function( name, extension, input.path = "input" ){
  file_name = paste(name, extension, sep = ".")
  file_path = file.path(dataset.path, file_name)
  if(!file.exists(file_path)) {
    Syn.input <- synStore(Syn.input)
    folder = Syn.input$properties$id
    
    syn.files = synChildren(folder)
    mask = names(syn.files) == file_name
    if(sum(mask) == 0) {
      warning(sprintf("Did not find file %s", file_name))
      return()
    } else {
      syn.id = as.character(syn.files[which(mask == T)[[1]]])
      entity <- synGet(syn.id, downloadLocation = dataset.path)  	
    }
  }
  if(extension == "csv") {
    obj = readr::read_csv(file_path, col_types=readr::cols())  
  } else if(extension == "tsv" | extension == "tab" | extension == "txt") {
    obj = readr::read_tsv(file_path, col_types=readr::cols())  
  } else if(extension == "rds" | extension == "RDS") {
    obj = readr::read_rds(file_path)
  }
  
  return(obj)
}

queryChEA3 <- function(genes, url = "https://maayanlab.cloud/chea3/api/enrich/") {
  library(httr)
  library(jsonlite)
  
  encode = "json"
  payload = list(query_name = "myQuery", gene_set = genes)
  
  #POST to ChEA3 server
  response = POST(url = url, body = payload, encode = encode)
  json = content(response, "text")
  
  #results as list of R dataframes
  results = fromJSON(json)
}
