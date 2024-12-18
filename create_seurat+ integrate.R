# Load the required libraries
library(Seurat)
library(Matrix)

# Define the parent directory containing the subfolders
parent_dir <- '/Users/mahimasrivastava/Downloads/cancer mastectomy 2'  # Replace with your path

# List all subfolders (samples)
sample_dirs <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)

# Initialize a list to store Seurat objects
seurat_objects <- list()

# Loop over each sample directory
for (sample_dir in sample_dirs) {
  # Read the matrix, barcodes, and features
  matrix <- readMM(file = file.path(sample_dir, "matrix.mtx.gz"))
  barcodes <- readLines(file.path(sample_dir, "barcodes.tsv.gz"))
  features <- read.table(file.path(sample_dir, "features.tsv.gz"), header = FALSE, sep = "\t")
  
  # Assign row names (genes) and column names (cells) to the matrix
  rownames(matrix) <- features$V1  # Gene names/IDs
  colnames(matrix) <- barcodes     # Cell barcodes
  
  # Create a Seurat object for this sample
  seurat_obj <- CreateSeuratObject(counts = matrix)
  
  # Add metadata (optional: you can modify or add sample-specific metadata)
  sample_name <- basename(sample_dir)  # Use folder name as sample name
  seurat_obj$sample <- sample_name
  
  # Store the Seurat object in the list
  seurat_objects[[sample_name]] <- seurat_obj
}

# Initialize the first Seurat object with its cell identifier
combined_seurat <- seurat_objects[[1]]
combined_seurat <- RenameCells(combined_seurat, add.cell.id = names(seurat_objects)[1])

# Loop over the remaining objects and merge them one by one
for (i in 2:length(seurat_objects)) {
  sample_name <- names(seurat_objects)[i]  # Use the sample name for the current object
  seurat_objects[[i]] <- RenameCells(seurat_objects[[i]], add.cell.id = sample_name)  # Add sample ID to the cells
  combined_seurat <- merge(combined_seurat, y = seurat_objects[[i]])
}

# Perform PCA
combined_seurat <- RunPCA(combined_seurat)

# Normalize and find variable features for each Seurat object individually
for (i in 1:length(seurat_objects)) {
  seurat_objects[[i]] <- NormalizeData(seurat_objects[[i]])
  seurat_objects[[i]] <- FindVariableFeatures(seurat_objects[[i]])
}

if (!requireNamespace("harmony", quietly = TRUE)) {
  devtools::install_github("immunogenomics/harmony")
}










# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)

# Define the parent directory
parent_dir <- '/Users/mahimasrivastava/Downloads/cancer mastectomy 2'

# List all subfolders (samples)
sample_dirs <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)

# Initialize a list to store Seurat objects
seurat_objects <- list()

# Loop over each sample directory
for (sample_dir in sample_dirs) {
  # Read matrix, barcodes, and features
  matrix <- readMM(file = file.path(sample_dir, "matrix.mtx.gz"))
  barcodes <- readLines(file.path(sample_dir, "barcodes.tsv.gz"))
  features <- read.table(file.path(sample_dir, "features.tsv.gz"), header = FALSE, sep = "\t")
  
  # Assign row and column names to the matrix
  rownames(matrix) <- features$V1
  colnames(matrix) <- barcodes
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = matrix)
  
  # Normalize and find variable features
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 5000)
  seurat_obj <- ScaleData(seurat_obj)
  
  # Add sample metadata
  sample_name <- basename(sample_dir)
  seurat_obj$sample <- sample_name
  
  # Store in the list
  seurat_objects[[sample_name]] <- seurat_obj
}

# Integration using canonical correlation analysis (CCA)
anchors <- FindIntegrationAnchors(object.list = seurat_objects, dims = 1:20, k.filter = 30, 
                                  anchor.features = 3000, k.score = 30)
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:20)

# Post-integration processing
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data)
integrated_data <- RunUMAP(integrated_data, dims = 1:20)
integrated_data <- FindNeighbors(integrated_data, dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 0.2)

# Visualization
DimPlot(integrated_data, reduction = "umap", group.by = "seurat_clusters")

# Save results
saveRDS(integrated_data, file = "integrated_data.rds")
























