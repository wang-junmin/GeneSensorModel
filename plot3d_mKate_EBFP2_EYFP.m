%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <- Author: Junmin Wang -> %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('/Users/junminwang/Documents/academics/PhDresearch/organoid_project/Fabio/codes/')

best_params = readtable('best_params_min_err_tbl.csv');
data = readtable('yfp_as_fn_of_binned_mKate_BFP_tidied_20x20.csv');

gRNA_level_lst = [1 2 3 4 8 12 16];
uORF_level_lst = [1 2];
for i=1:7
    gRNA_level = gRNA_level_lst(i);
    for j=1:2
        uORF_level = uORF_level_lst(j);
        tmp_data = data(data.gRNA == gRNA_level & ...
                data.uORF == uORF_level, :);
        EBFP2_level_lst = unique(tmp_data.EBFP2);
        for m = 1:length(EBFP2_level_lst)
            EBFP2_level = EBFP2_level_lst(m);
            tmp_data_2 = tmp_data(tmp_data.EBFP2 == EBFP2_level, :);
            p1 = plot3(tmp_data_2.mKate, ...
                tmp_data_2.EBFP2, ...
                tmp_data_2.EYFP, ...
                'o', 'Color', ...
                hsv2rgb([m/(length(EBFP2_level_lst)+1) 1 1]), ...
                'MarkerFaceColor', ...
                hsv2rgb([m/(length(EBFP2_level_lst)+1) 1 1]), ...
                'MarkerSize', 3);
            hold on
        end
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        set(gca, 'ZScale', 'log')  
        set(gca,'XMinorTick', 'on', ...
            'YMinorTick', 'on', ...
            'ZMinorTick', 'on')
        xlim([10^3 10^6]);
        ylim([10^2 10^6.7]);
        zlim([10^4.5 10^8.5]);
        xlabel('mKate [MEFL]', 'FontSize', 20);
        ylabel('EBFP2 [MEFL]', 'FontSize', 20);
        zlabel('EYFP [MEFL]', 'FontSize', 20);
        set(gca, 'fontsize', 15);
        [X,Y] = meshgrid(2:0.5:6, 2:0.5:6);
        Z = best_params.k0 + (best_params.c * gRNA_level / (gRNA_level + best_params.m) + best_params.d * uORF_level) ...
            ./ (1 + exp(best_params.a0 + best_params.a1 * X + best_params.a2 * Y + best_params.a12 .* X .* Y ...
            + best_params.a3 * gRNA_level + best_params.a4 * uORF_level));
        s = surf(10.^X, 10.^Y, 10.^Z); %
        set(s,'FaceColor', [0.5843 0.8157 0.9882], ...
            'FaceAlpha', 0.1);
        s.EdgeColor = 'none';
        savefig(sprintf('mKate_EBFP2_EYFP_gRNA_%d_uORF_%d.fig', ...
            gRNA_level, uORF_level));
        close
    end
end

