import argparse
import os
import pickle
import sys
import textwrap

import celltypist
import scanpy as sc


def train(input_path,
          labels,
          genes,
          model_name,
          output_path,
          transpose_input,
          normalize,
          raw_counts_layer,
          gene_symbols_field,
          epochs,
          max_iter,
          C,
          n_jobs,
          use_SGD,
          alpha,
          use_GPU,
          mini_batch,
          batch_number,
          batch_size,
          balance_cell_type,
          feature_selection,
          top_genes):
    """
    Method used to train a CellTypist model and output the model to a desired location. 
    Please refer to the CLI help to gain a better understanding of the required and optional
    parameters.
    """
    if input_path.endswith('.h5ad'):
        print("Detected an AnnData object as input.")
        data = sc.read_h5ad(input_path)
    else:
        print("Detected .mtx file as input.")
        data = input_path
    if input_path.endswith('.mtx') and genes == None:
        raise ValueError("Missing a gene file for the provided data.")
    if raw_counts_layer and input_path.endswith('.h5ad') and raw_counts_layer in data.layers.keys() and raw_counts_layer!= 'X':         
        print(f"Using raw counts layer: {raw_counts_layer}")
        data.X = data.layers[raw_counts_layer]
    elif raw_counts_layer and input_path.endswith('.h5ad'):
        raise ValueError(f"Raw counts layer {raw_counts_layer} should be either different to 'X' or an existing layer in AnnData")
    elif raw_counts_layer:
        raise ValueError(f"Raw counts layer {raw_counts_layer} provided but the data format provided is .mtx. Please provide an AnnData object.")
    if gene_symbols_field and input_path.endswith('.h5ad') and gene_symbols_field in data.var.columns:
        data.var['old_index'] = data.var.index
        data.var.set_index(gene_symbols_field, inplace=True)
    if normalize:
        sc.pp.normalize_total(data, target_sum=1e4)
        sc.pp.log1p(data)
        
    model = celltypist.train(data,
                             labels=labels,
                             genes=genes,
                             transpose_input=transpose_input,
                             C=C,
                             epochs=epochs,
                             max_iter=max_iter,
                             n_jobs=n_jobs,
                             use_SGD=use_SGD,
                             alpha=alpha,
                             use_GPU=use_GPU,
                             mini_batch=mini_batch,
                             batch_number=batch_number,
                             batch_size=batch_size,
                             balance_cell_type=balance_cell_type,
                             feature_selection=feature_selection,
                             top_genes=top_genes)
    
    model_save_location = output_path + '/' + model_name + '.pkl'
    
    return pickle.dump(model, open(model_save_location, 'wb'))

def predict(input_path, 
            model, 
            output_path,
            gene_file,
            cell_file,
            raw_counts_layer,
            gene_symbols_field,
            normalize, 
            transpose_input,
            mode, 
            p_thres,
            majority_voting,
            over_clustering,
            use_GPU):
    """
    Method used to obtain cell type prediction on a new dataset using an existing CellTypist
    model.
    Please refer to the CLI help to gain a better understanding of the required and optional
    parameters.
    """
    if input_path.endswith('.h5ad'):
        print("Detected an AnnData object as input.")
        data = sc.read_h5ad(input_path)
    else:
        print("Detected .mtx file as input.")
        data = input_path
    if model not in celltypist.models.get_all_models():
        with open(model, 'rb') as pickle_load:
            model = pickle.load(pickle_load)
    if input_path.endswith('.mtx') and gene_file == None:
        raise ValueError("Missing a gene file for the provided .mtx data.")
    if input_path.endswith('.mtx') and cell_file == None:
        raise ValueError("Missing a cell file for the provided .mtx data.")
    if raw_counts_layer and input_path.endswith('.h5ad') and raw_counts_layer in data.layers.keys() and raw_counts_layer!= 'X':         
        print(f"Using raw counts layer: {raw_counts_layer}")
        data.X = data.layers[raw_counts_layer]
    elif raw_counts_layer and input_path.endswith('.h5ad'):
        raise ValueError(f"Raw counts layer {raw_counts_layer} should be either different to 'X' or an existing layer in AnnData")
    elif raw_counts_layer:
        raise ValueError(f"Raw counts layer {raw_counts_layer} provided but the data format is .mtx. Please provide an AnnData object.")
    if gene_symbols_field and input_path.endswith('.h5ad') and gene_symbols_field in data.var.columns:
        data.var['old_index'] = data.var.index
        data.var.set_index(gene_symbols_field, inplace=True)
    if normalize:
        sc.pp.normalize_total(data, target_sum=1e4)
        sc.pp.log1p(data)
    if mode:
        mode = mode.replace('_', ' ')
    if majority_voting:
        if over_clustering:
            over_clustering = over_clustering
        else:
            over_clustering = None
            
    predictions = celltypist.annotate(data,
                                      model,
                                      gene_file=gene_file,
                                      cell_file=cell_file,
                                      mode=mode, 
                                      p_thres=p_thres, 
                                      transpose_input=transpose_input,
                                      majority_voting=majority_voting,
                                      over_clustering=over_clustering,
                                      use_GPU=use_GPU)
    
    return predictions.predicted_labels.to_csv(output_path)

def main():
    parser = argparse.ArgumentParser(add_help=False,
                                     description=textwrap.dedent('''
                                     CellTypist for automatic cell type annotation.
                                     ----------------------------------------------
                                     Welcome to the ODS CLI tool for cell type 
                                     annotation using CellTypist. 
                                     The CLI provides the user with the opportunity 
                                     to train a new model or predict using an 
                                     existing model. 
                                     Either of the methods can be activated by 
                                     providing the appropriate action variable.
                                     For details of required/optional parameters
                                     for each method, please refer to the 
                                     appropriate help function (-h).
                                     '''),
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    choices_to = {"train": "⚖️  Train a CellTypist model using a reference dataset.",
                 "predict": "🖋️  Use an existing CellTypist model to predict on a new dataset."}
    
    parser.add_argument("--action",
                       type=str,
                       help='\n'.join("{}: {}".format(key, value) for key, value in choices_to.items()),
                       choices=choices_to)
    
    parser.add_argument('--help', action='store_true', help="List available options.", required=False)
    
    args, sub_args = parser.parse_known_args()
    if args.help:
        if args.action is None: 
            print(parser.format_help())
            sys.exit(1)
        sub_args.append('--help')
        
    action = "predict" if args.action is None else args.action

    parser = argparse.ArgumentParser(prog="%s %s" % (os.path.basename(sys.argv[0]), action))
    
    
    if action == "predict":
        parser.add_argument("--input_path", 
                            type=str, 
                            help="Path to the input file. Can be an AnnData object or a .mtx file.", 
                            required=True)

        parser.add_argument("--output_path", 
                            type=str, 
                            help="Path to output file.", 
                            required=True)

        parser.add_argument("--model", 
                            type=str, 
                            help="Path to model file, stored in .pkl format.", 
                            required=True)
        
        parser.add_argument("--gene_file",
                           type=str,
                           help="Path to the file which stores each gene per line corresponding to the genes used in the provided .mtx file.",
                           required=False,
                           default=None)
        
        
        parser.add_argument("--cell_file",
                           type=str,
                           help="Path to the file which stores each cell per line corresponding to the cells used in the provided .mtx file.",
                           required=False,
                           default=None)

        parser.add_argument("--raw_counts_layer", 
                            type=str, 
                            help="The name of the layer that stores the raw counts. Uses default matrix if not present.",
                            required=False, 
                            default=None)

        parser.add_argument("--gene_symbols_field", 
                            type=str, 
                            help="The field in AnnData where the gene symbols are stored, if not in index.", 
                            required=False, 
                            default=None)

        parser.add_argument("--normalize", 
                            action='store_true', 
                            help="If raw counts are provided in the AnnData object, they need to be normalized.", 
                            required=False)

        parser.add_argument("--transpose_input", 
                            action='store_true', 
                            help="If the provided matrix is in the gene-by-cell format, please transpose the input to cell-by-gene format", 
                            required=False)

        parser.add_argument("--mode", 
                            type=str, 
                            help="Mode for the prediction is to choose the cell type, defualts to best match.", 
                            required=False, 
                            default='best match')

        parser.add_argument("--p_thres", 
                            type=float, 
                            help="Probability threshold for assigning a cell type in a multiclass problem, defualts to 0.5.",
                            default=0.5, 
                            required=False)

        parser.add_argument("--majority_voting", 
                            action='store_true', 
                            help="Refine the predicted labels by running the majority voting classifier after over-clustering.",
                            required=False)

        parser.add_argument("--over_clustering", 
                            type=str, 
                            help="If majority voting is set to True, specify the type of over clustring that is to be perfomend. This can be specified in the AnnData or an input file specifying the over-clustering per cell. If not present, then the default heuristic over-clustring based on input data will be used.", 
                            required=False, 
                            default=None)

        parser.add_argument("--use_GPU", 
                            action='store_true', 
                            help="Whether to use GPU for over clustering on the basis of `rapids-singlecell`.", 
                            required=False)


        
        args=parser.parse_args(sub_args)
        
        predict(args.input_path, 
                args.model, 
                args.output_path,
                args.gene_file,
                args.cell_file,
                args.raw_counts_layer, 
                args.gene_symbols_field, 
                args.normalize, 
                args.transpose_input, 
                args.mode, 
                args.p_thres,
                args.majority_voting, 
                args.over_clustering, 
                args.use_GPU)
    else:
        parser.add_argument("--input_path", 
                            type=str, 
                            help="Path to the input file. Can be an AnnData object or a .mtx file.", 
                            required=True)
        
        parser.add_argument("--labels",
                           type=str,
                           help="Path to the file that stores the per cell types or the layer in the AnnData object where the cell types are held.",
                           required=True)
        
        parser.add_argument("--genes",
                           type=str,
                           help="Path to the file containing one gene per line corresponding to the genes in X, required for .mtx data.",
                           required=False,
                           default=None)
        
        parser.add_argument("--model_name",
                           type=str,
                           help="The name of the trained model, used for saving.",
                           required=True)
        
        parser.add_argument("--output_path",
                           type=str,
                           help="The location to where the trained model should be saved.",
                           required=True)
        
        parser.add_argument("--transpose_input", 
                            action='store_true', 
                            help="If the provided matrix is in the gene-by-cell format, please transpose the input to cell-by-gene format.", 
                            required=False)
        
        parser.add_argument("--normalize", 
                            action='store_true', 
                            help="If raw counts are provided in the AnnData object, they need to be normalized.", 
                            required=False)
        
        parser.add_argument("--gene_symbols_field", 
                            type=str, 
                            help="The field in AnnData where the gene symbols are stored, if not in index.", 
                            required=False, 
                            default=None)
        
        parser.add_argument("--raw_counts_layer", 
                            type=str, 
                            help="The name of the layer that stores the raw counts. Uses default matrix if not present",
                            required=False, 
                            default=None)
        
        parser.add_argument("--epochs", 
                            type=int, 
                            help="The number of epochs for which the model needs to be trained.",
                            required=True, 
                            default=10)
        
        parser.add_argument("--solver", 
                            type=str, 
                            help="Algorithm to use in the optimization problem for traditional logistic classifier. Default is based on the size of the input data.",
                            required=False)
        
        parser.add_argument("--max_iter", 
                            type=int, 
                            help="Maximum number of iterations before reaching the minimum of the cost function.",
                            required=True, 
                            default=100)
        
        parser.add_argument("--C", 
                            type=float, 
                            help="Inverse of L2 regularization strength for traditional logistic classifier.", 
                            required=False, 
                            default=1.0)
        
        parser.add_argument("--n_jobs", 
                            type=int, 
                            help="Number of CPUs used.", 
                            required=False, 
                            default=1)
        
        parser.add_argument("--use_SGD", 
                            type=bool, 
                            help="Whether to implement SGD learning for the logistic classifier.", 
                            required=False, 
                            default=False)
        
        parser.add_argument("--alpha", 
                            type=float, 
                            help="L2 regularization strength for SGD logistic classifier.", 
                            required=False, 
                            default=0.0001)
        
        parser.add_argument("--use_GPU", 
                            type=bool, 
                            help="Whether to use GPU for logistic classifier.", 
                            required=False, 
                            default=False)
        
        parser.add_argument("--mini_batch", 
                            type=bool, 
                            help="Whether to implement mini-batch training for the SGD logistic classifier.", 
                            required=False, 
                            default=False)
        
        parser.add_argument("--batch_number", 
                            type=int, 
                            help="The number of batches used for training in each epoch. Each batch contains batch_size cells.", 
                            required=False, 
                            default=100)
        
        parser.add_argument("--batch_size", 
                            type=int, 
                            help="The number of cells within each batch.", 
                            required=False, 
                            default=1000)
        
        parser.add_argument("--balance_cell_type", 
                            type=bool, 
                            help="Whether to balance the cell type frequencies in mini-batches during each epoch.", 
                            required=False, 
                            default=False)
        
        parser.add_argument("--feature_selection", 
                            type=bool, 
                            help="Whether to perform two-pass data training where the first round is used for selecting important features/genes using SGD learning.", 
                            required=False, 
                            default=False)
        
        parser.add_argument("--top_genes", 
                            type=int, 
                            help="The number of top genes selected from each class/cell-type based on their absolute regression coefficients.", 
                            required=False, 
                            default=300)
        
        args = parser.parse_args(sub_args)
        
        train(args.input_path,
              args.labels,
              args.genes,
              args.model_name,
              args.output_path,
              args.transpose_input,
              args.normalize,
              args.raw_counts_layer,
              args.gene_symbols_field,
              args.epochs,
              args.max_iter,
              args.C,
              args.n_jobs,
              args.use_SGD,
              args.alpha,
              args.use_GPU,
              args.mini_batch,
              args.batch_number,
              args.batch_size,
              args.balance_cell_type,
              args.feature_selection,
              args.top_genes)

if __name__ == "__main__":
    main()

    
