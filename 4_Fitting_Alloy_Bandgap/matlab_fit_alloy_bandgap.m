clc, clear, close all
warning('off')

[E221, dos221] = load_dos('dos-221-lro0-sro0-k444.dat', 64);    % filename, atom number
[E444, dos444] = load_dos('dos-442-lro0-sro0-k111.dat', 512);

figure(1)
set(gcf, 'position', [100 100 600 670])
axes('position', [0.16 0.6 0.7 0.3])
hold on, box on
[E_vbm_221, E_cbm_221] = get_fit_dos(E221, dos221, 1, -0.3, 0.5, -1.4); % cbm_start, vbm_start, cbm_range, vbm_range
set(gca, 'xtick', -5:1:5, 'ytick', 0:1000:2000)

axes('position', [0.16 0.21 0.7 0.3])
hold on, box on
[E_vbm_444, E_cbm_444] = get_fit_dos(E444, dos444, 0.8, -0.6, 0.8, -1.5);
xlabel('Energy (eV)')
set(gca, 'xtick', -5:1:5, 'ytick', 0:1000:2000)




function [loc_E_0_vbm, loc_E_0_cbm] = get_fit_dos(E,  dos, st_cbm, st_vbm, cbm_step, vbm_step)

%%% plot DOS
plot(E, dos, '-', 'linewidth', 1, 'color', [0 152 255]/255)

%%% fit CBM
[fit_E, fit_dos] = get_fit_data(E, dos, st_cbm, cbm_step, 'cbm');
[Ei, cur_doii] = fit_equation(fit_E, fit_dos);
[~, doi_0_id] = min(abs(cur_doii));
loc_E_0_cbm = Ei(doi_0_id);

plot(fit_E, fit_dos, 'ok', 'linewidth', 1)
plot(Ei, cur_doii, 'r-', 'linewidth', 1)


%%% fit VBM
[fit_E, fit_dos] = get_fit_data(E, dos, st_vbm, vbm_step, 'vbm');
[Ei, cur_doii] = fit_equation(-fit_E, fit_dos);
[~, doi_0_id] = min(abs(cur_doii ));
loc_E_0_vbm = -Ei(doi_0_id);
plot(fit_E, fit_dos, 'ok', 'linewidth', 1)
plot(-Ei, cur_doii, 'r-', 'linewidth', 1)
xlim([-5 5])
ylim([0 2000])
set(gca, 'fontname', 'arial', 'fontsize', 12, 'linewidth', 0.8, ...
    'xcolor', [0 0 0], 'ycolor', [0 0 0])

%%% bandgap
bandgap = loc_E_0_cbm - loc_E_0_vbm

end

function [Ei, dosi] = get_fit_data(E, dos, st, step, loc)
if strcmp(loc, 'vbm')
    fit_E_cur = linspace(st+step, st, 21);
else 
    if strcmp(loc, 'cbm')
        fit_E_cur = linspace(st, st+step, 21);
    end
end

fit_Eid = zeros(size(fit_E_cur));
for j = 1:length(fit_E_cur)
    [~, nearest_id] = min(abs(E-fit_E_cur(j)));
    fit_Eid(j) = nearest_id;
end
Ei = E(fit_Eid)';
dosi = dos(fit_Eid);
end

function [xi, yi] = fit_equation(x, y)
f = fittype('A.*(x+1./E.*(5.*x.^2)+1./(E.^2).*(8.*x.^3)+1./(E.^3).*(4.*x.^4))+C', 'independent', 'x', 'coefficients', {'A', 'E', 'C'});
cfun = fit(x', y', f, 'lower', [-1 -1 -1]*100000, 'upper', [1 1 1]*1000000);
xi = linspace(min(x)-0.5, max(x)+0.5, 1e3);
yi = cfun(xi);
end

function [E, dos] = load_dos(file_name, Natom)
data = importdata(file_name);
data = data.data;
E = data(:,1);
dos = (data(:,2)'*64/Natom).^2;
end