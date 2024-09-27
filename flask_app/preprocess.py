import os
import dandelion as ddl
import scanpy as sc

def preprocess(upload_folder, samples, data_uploaded, species = 'human'):
    sample_paths = [os.path.join(upload_folder, sample) for sample in samples]
    
    if data_uploaded == 'Both':
    
        ddl.pp.format_fastas(sample_paths, prefix=samples)
        ddl.pp.reannotate_genes(sample_paths, org = species)
        ddl.pp.assign_isotypes(sample_paths, plot = False, org = species)
        
        adata_list = []
        for sample in sample_paths:
            # Find the .h5 file in the sample folder
            h5_file = None
            for file in os.listdir(sample):
                if file.endswith('.h5'):
                    h5_file = file
                    break
            
            if h5_file is None:
                print(f"No .h5 file found in {sample}")
                continue
            
            h5_file_path = os.path.join(sample, h5_file)
            
            adata = sc.read_10x_h5(
                h5_file_path,
                gex_only=True,
            )
            adata.obs["sampleid"] = sample
            # rename cells to sample id + barcode, and cleave the trailing -#
            adata.obs_names = [
                str(sample) + "_" + str(j).split("-")[0] for j in adata.obs_names
            ]
            adata.var_names_make_unique()
            adata_list.append(adata)
        # no need for index_unique as we already made barcodes unique by prepending the sample ID
        adata = adata_list[0].concatenate(adata_list[1:], index_unique=None)
        
        vdj_list = []
        for sample in sample_paths:
            tsv_file = None
            for file in os.listdir(os.path.join(sample, 'dandelion')):
                if file.endswith('.tsv'):
                    tsv_file = 'dandelion/' + file
                    break
            
            if tsv_file is None:
                print(f"No .tsv file found in {sample}")
                continue
            
            tsv_file_path = os.path.join(sample, tsv_file)
            
            vdj = ddl.read_10x_airr(
                tsv_file_path
            )
            # the dandelion output already has the sample ID prepended at the start of each contig
            vdj_list.append(vdj)
            
        vdj = ddl.concat(vdj_list)
        
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver="arpack")
        sc.pp.neighbors(adata, n_pcs=20)
        sc.tl.umap(adata)
        
        sc.tl.leiden(adata, resolution=0.5)
        
        vdj, adata = ddl.pp.check_contigs(vdj, adata)
        ddl.tl.find_clones(vdj)
        ddl.tl.generate_network(vdj)
        ddl.tl.clone_size(vdj)
        
        
        adata.write(os.path.join(upload_folder, "processed_adata.h5ad"))
        vdj.write(os.path.join(upload_folder, "processed_vdj.h5ddl"))
        
        return os.path.join(upload_folder, "processed_adata.h5ad"), os.path.join(upload_folder, "processed_vdj.h5ddl")
    
    else:
        
        ddl.pp.format_fastas(sample_paths, prefix=samples)
        ddl.pp.reannotate_genes(sample_paths, org = species)
        ddl.pp.assign_isotypes(sample_paths, plot = False, org = species)
        
        vdj_list = []
        for sample in sample_paths:
            tsv_file = None
            for file in os.listdir(os.path.join(sample, 'dandelion')):
                if file.endswith('.tsv'):
                    tsv_file = 'dandelion/' + file
                    break
            
            if tsv_file is None:
                print(f"No .tsv file found in {sample}")
                continue
            
            tsv_file_path = os.path.join(sample, tsv_file)
            
            vdj = ddl.read_10x_airr(
                tsv_file_path
            )
            # the dandelion output already has the sample ID prepended at the start of each contig
        vdj_list.append(vdj)
        vdj = ddl.concat(vdj_list)
        vdj = ddl.pp.check_contigs(vdj)
        ddl.tl.find_clones(vdj)
        ddl.tl.generate_network(vdj)
        ddl.tl.clone_size(vdj)
            
        vdj.write(os.path.join(upload_folder, "processed_vdj.h5ddl"))
        return os.path.join(upload_folder, "processed_vdj.h5ddl")

    
if __name__ == "__main__":
    preprocess('/Users/nick/uploads/test', samples = ["sc5p_v2_hs_PBMC_1k","sc5p_v2_hs_PBMC_10k"])


