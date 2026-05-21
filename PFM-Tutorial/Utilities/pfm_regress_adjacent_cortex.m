function [Output] = pfm_regress_adjacent_cortex(Input,D,Distance)

% count the number of cortical vertices (should be 59412);
nCorticalVertices = nnz(Input.brainstructure==1) + nnz(Input.brainstructure==2);

% load distance matrix;
if ischar(D)
    D = smartload(D); 
end

% index of subcortical voxels
SubcortVoxels = (nCorticalVertices+1):size(D,1);

% preallocate;
Output = Input;

% sweep all
% subcortical voxels
% nearby gray matter;
% extract nearby gm signals;
% average; if needed
% remove (possible) contamination of nearby cortical signals via linear regression

nb_mask = sparse(D(1:nCorticalVertices,SubcortVoxels) <= Distance);
wh_nb_sc = any(nb_mask);
nb_mask = nb_mask(:, wh_nb_sc);

nb_gm_ts = Input.data(1:nCorticalVertices,:).' * nb_mask ./ sum(nb_mask);
nb_gm_ts = nb_gm_ts - mean(nb_gm_ts);
nb_sc_ts = Input.data(nCorticalVertices+find(wh_nb_sc),:).';
nb_sc_ts = nb_sc_ts - mean(nb_sc_ts);
nb_gm_b = sum(nb_gm_ts .* nb_sc_ts) ./ sum(nb_gm_ts.^2);
nb_sc_resid = nb_sc_ts - nb_gm_b .* nb_gm_ts;

Output.data(nCorticalVertices+find(wh_nb_sc),:) = nb_sc_resid.';

end

% subfunctions
function out = smartload(matfile)
out = load(matfile);
names = fieldnames(out);
out = eval(['out.' names{1}]);
end