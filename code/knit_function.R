#################
# KNIT FUNCTION #
#################
custom_knit <- function(inputFile, out_dir, ...) {
     rmarkdown::render(
      inputFile,
      output_file = paste0(
        xfun::sans_ext(inputFile), '_', Sys.Date(),'.html'),
      output_dir = paste0("../lab_book/",out_dir,"/"),
      output_format = rmarkdown::html_document(
        dev = c("jpeg", "tiff", "pdf"),
        dpi = 300,
        echo=TRUE,
        #keep_md = TRUE,
        toc = FALSE,                             
        toc_float = FALSE,
        code_folding = "show",
        code_download = TRUE)
      )#;
    # rmarkdown::render(
    #   inputFile,
    #   output_format = rmarkdown::pdf_document()
    #   )
}

##################
# KNIT RMARKDOWN #
#################
#### html ####
knit_gfm_with_job <- function() {
  rstudioapi::verifyAvailable()
  
  job_file <- tempfile(fileext = ".R")
  
  active_doc_ctx <- rstudioapi::getSourceEditorContext()
  rmd_path <- active_doc_ctx$path
  
  if (identical(rmd_path, "")) {
    rstudioapi::showDialog(
      "Cannot Render Unsaved R Markdown Document",
      "Please save the current document before rendering."
    )
    return(invisible())
  }
  
  rstudioapi::documentSave(active_doc_ctx$id)
  
  rmd_path <- normalizePath(rmd_path, mustWork = TRUE)
  out_file <- paste0(basename(xfun::sans_ext(rmd_path)), '.md')
  
  cat(
    'res <- rmarkdown::render("', rmd_path, '",',
    '"', out_file, '",',
    ' output_format = rmarkdown::github_document(html_preview = FALSE, dev = "jpeg")', ')\n',
    'unlink("', job_file, '")\n',
    'rstudioapi::viewer(res)',
    sep = "",
    file = job_file
  )
  
  rstudioapi::jobRunScript(
    path = job_file,
    name = basename(rmd_path),
    workingDir = dirname(rmd_path),#out_dir,
    importEnv = FALSE
  )
} 

#params = list(print=TRUE, setup=TRUE),

# https://stackoverflow.com/questions/66620582/customize-the-knit-button-with-source
########################################
# CONVERT .MD/.TEX TO MULTIPLE FORMATS #
########################################

# knit: (function(inputFile, encoding) {
#      rmarkdown::render(
#       inputFile,
#       encoding = encoding,
#       params = list(print=TRUE, setup=TRUE),
#       output_format = 'md_document'
#       );
#       purrr::walk(c("pdf","docx", "html"), 
#                 ~rmarkdown::pandoc_convert(
#                     "Figure1.tex", 
#                     output=paste0("Figure1.", .x)))                      
#                                 })

##############################
# CONDITIONAL GLOBAL OPTIONS #
##############################

# if(knitr::is_latex_output()){
#   knitr::opts_chunk$set(
#     echo=FALSE,
#     dev=c('pdf', 'png'),
#     dpi = 300)
# }
# 
# if(isFALSE(knitr::is_latex_output())){
#   knitr::opts_chunk$set(
#     echo=TRUE,
#     toc = FALSE,
#     toc_float = FALSE,
#     code_folding = "show",
#     code_download = TRUE)
# }

#######################
# TEST YAML FORMATING #
#######################
# cat(yaml::as.yaml(list(
#   title = 'A Wonderful Day',
#   dev = c('png', 'pdf'),
#   dpi = 300
# )))
##################
# YAML HEADER #
##################
# ---
# title: "Figure 1.Luminal Study Groups"
# geometry: "left=2cm,right=2cm,top=2cm,bottom=2cm"
# output: 
#   md_document:
#     preserve_yaml: true
#   pdf_document:
#     dev:
#     - png
#     - pdf
#     fig_caption: yes
#   html_document:
#     code_download: true
#     code_folding: show
#     dpi: 300.0
#     toc: yes
#     toc_float:
#       collapsed: no
# header-includes: 
# - \usepackage{float}
# editor_options: 
#   chunk_output_type: console
# params:
#   setup: TRUE
#   print: FALSE