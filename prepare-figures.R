# This script just reorganizes all figures for presentation in one folder, with meaningful names for the manuscript


# Preprint ----------------------------------------------------------------

figure_correspondence_table = read_excel("./other-data/Manuscript Figure Correspondence.xlsx", sheet = 1)

output_path = paste0("./output/", results_path)

paper_fig_path = paste0("./output/", results_path, "/_manuscript figures and tables/preprint/")
suppressWarnings({dir.create(paper_fig_path)})
suppressWarnings({dir.create(paste0(paper_fig_path, "/tables"))})


message(paste0("Saving figures and tables to ", results_path))

pwalk(figure_correspondence_table, function (code_generated_filename, manuscript_name, type, ...) {
  fnf = paste0(output_path, code_generated_filename)
  fnt = paste0(paper_fig_path, str_to_title(type), "_", manuscript_name, str_extract(basename(fnf), "\\..+?$"))
  
  if(type == "table"){
    fnt = paste0(paper_fig_path, "/tables/", str_to_title(type), " ", manuscript_name, str_extract(basename(fnf), "\\..+?$"))
  } 

  if (!is.na(code_generated_filename)) {
    if (file.exists(fnf)) {
    file.copy(from = fnf, to = fnt, overwrite = T)
    message(paste0(basename(fnt), " <-- ", fnf))
    } else {
      message(basename(fnt), paste0(" <-- ERROR - file not found: ", fnf))
    }
  } else {
    message(paste0(basename(fnt), " <-- SKIP"))
  }
})


# Nature ------------------------------------------------------------------

figure_correspondence_table = read_excel("./other-data/Manuscript Figure Correspondence.xlsx", sheet = 2)

output_path = paste0("./output/", results_path)

paper_fig_path = paste0("./output/", results_path, "/_manuscript figures and tables/nature/")
suppressWarnings({dir.create(paper_fig_path)})
suppressWarnings({dir.create(paste0(paper_fig_path, "/tables"))})

message(paste0("Saving figures and tables to ", results_path))

pwalk(figure_correspondence_table, function (code_generated_filename, manuscript_name, type, ...) {
  fnf = paste0(output_path, code_generated_filename)
  fnt = paste0(paper_fig_path, manuscript_name, str_extract(basename(fnf), "\\..+?$"))
  
  if(type == "table"){
    fnt = paste0(paper_fig_path, "/tables/", manuscript_name, str_extract(basename(fnf), "\\..+?$"))
  } 
  
  if (!is.na(code_generated_filename)) {
    if (file.exists(fnf)) {
      file.copy(from = fnf, to = fnt, overwrite = T)
      message(paste0(basename(fnt), " <-- ", fnf))
    } else {
      message(basename(fnt), paste0(" <-- ERROR - file not found: ", fnf))
    }
  } else {
    message(paste0(basename(fnt), " <-- SKIP"))
  }
})

# Deleting tables folder
file.remove(list.files(path = paste0(output_path, "/_manuscript figures and tables/tables"),
                       pattern = "Table[ _](S\\d+|[123])",
                       full.names = TRUE,
                       ignore.case = TRUE))

file.rename(paste0(output_path, "/_manuscript figures and tables/tables"), paste0(output_path, "/_manuscript figures and tables/other tables"))
