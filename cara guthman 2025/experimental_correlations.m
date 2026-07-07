function [ks_stat, ks_pval, corrTbl] = experimental_correlations(mpo, vmh, filename)
%inits
E = height(mpo);
comparisons = combinations(1:E, 1:E);
headers = {'rho', 'MPO_Exp', 'VMH_Exp', 'Type'};
vartypes = {'double', 'categorical', 'categorical', 'categorical'};
corrTbl = table('Size', [0 numel(headers)], 'VariableTypes', vartypes, 'VariableNames', headers);

%replace nans with noise at z = 0
mpo(isnan(mpo)) = eps; vmh(isnan(vmh)) = eps;

for comp = 1:height(comparisons)
    %initialize comparison
    e1 = comparisons.Var1(comp); e2 = comparisons.Var2(comp);

    %pull comparison and put into table
    tmp = corrcoef(mpo(e1, :), vmh(e2, :));
    if e1 == e2
        corrTbl = cat(1, corrTbl, {tmp(2), categorical(e1), categorical(e2), "matched"});
    else
        corrTbl = cat(1, corrTbl, {tmp(2), categorical(e1), categorical(e2), "mismatched"});
    end
end

%get data and stats
matched = corrTbl.rho(corrTbl.Type == "matched", :);
[fm, xm] = ecdf(matched);
mismatched = corrTbl.rho(corrTbl.Type == "mismatched", :);
[fmm, xmm] = ecdf(mismatched);
[~, ks_pval, ks_stat] = kstest2(matched, mismatched);

%plot
ksFig = figure("Position", [400 400 320 320], "Color", [1 1 1]);
hold on
plot(xmm, fmm, 'Color', [0 0 0])
plot(xm, fm, 'Color', [62 162 44] ./ 255)

%aesthetics
ksAx = gca;
setAx(ksAx);
axis square
xlim([-1 1])
xlabel("correlation coef")
ylim([0 1])
ylabel("cumulative prob")

%save and close
exportgraphics(ksFig, filename, "ContentType", "vector")
close all
end