clc, clear, close all
warning('off')

[E221, dos221] = load_dos('dos-221-lro0-sro0-k444.dat', 64);
[E222, dos222] = load_dos('dos-222-lro0-sro0-k331.dat', 128);
[E331, dos331] = load_dos('dos-331-lro0-sro0-k334.dat', 144);
[E332, dos332] = load_dos('dos-332-lro0-sro0-k221-1.dat', 288);
[E444, dos444] = load_dos('dos-442-lro0-sro0-k111.dat', 512);
[E554, dos554] = load_dos('dos-552-lro0-sro0-k111.dat', 800);

figure(1)
set(gcf, 'position', [100 100 600 670])
axes('position', [0.16 0.77 0.7 0.13])
hold on, box on
[E_vbm_221, E_cbm_221] = get_fit_dos(E221, dos221, 1, -0.3, 0.5, -1.4);
set(gca, 'xtick', -5:1:5, 'ytick', 0:1000:2000)

axes('position', [0.16 0.63 0.7 0.13])
hold on, box on
[E_vbm_222, E_cbm_222] = get_fit_dos(E222, dos222, 0.45, -0.65, 0.9, -1.3);
set(gca, 'xtick', -5:1:5, 'ytick', 0:1000:2000)

axes('position', [0.16 0.49 0.7 0.13])
hold on, box on
[E_vbm_331, E_cbm_331] = get_fit_dos(E331, dos331, 0.6, -0.5, 0.4, -1.5);
ylabel('DOS^2')
set(gca, 'xticklabel', {}, 'xtick', -5:1:5, 'ytick', 0:1000:2000)
% 
axes('position', [0.16 0.35 0.7 0.13])
hold on, box on
[E_vbm_332, E_cbm_332] = get_fit_dos(E332, dos332, 0.5, -0.9, 0.6, -1.5);
set(gca, 'xtick', -5:1:5, 'ytick', 0:1000:2000)

% 0.7, -0.7, 1.1, -1.5 : 2
axes('position', [0.16 0.21 0.7 0.13])
hold on, box on
[E_vbm_444, E_cbm_444] = get_fit_dos(E444, dos444, 0.8, -0.6, 0.8, -1.5);
xlabel('Energy (eV)')
set(gca, 'xtick', -5:1:5, 'ytick', 0:1000:2000, 'xticklabel', {})

% 
axes('position', [0.16 0.07 0.7 0.13])
hold on, box on
[E_vbm_444, E_cbm_444] = get_fit_dos(E554, dos554, 0.65, -0.6, 1.2, -1.5);
xlabel('Energy (eV)')
set(gca, 'xtick', -5:1:5, 'ytick', 0:1000:2000)





function [loc_E_0_vbm, loc_E_0_cbm] = get_fit_dos(E,  dos, st_cbm, st_vbm, cbm_step, vbm_step)
%%% cbm
fit_E = [];
fit_dos = [];

E_pos_id = E > st_cbm & E < 3;
dos_pos = dos(E_pos_id);
E_pos = E(E_pos_id);
[~, E_pos_id] = min(abs(E_pos-st_cbm));
% [~, E_pos_id] = min(abs(dos_pos-0.001));
E_start = E_pos(E_pos_id);
fit_E_cur = linspace(E_start, E_start+cbm_step, 21);
fit_Eid = zeros(size(fit_E_cur));
for j = 1:length(fit_E_cur)
    [~, nearest_id] = min(abs(E-fit_E_cur(j)));
    fit_Eid(j) = nearest_id;
end
fit_E = [fit_E; E(fit_Eid)'];
fit_dos = [fit_dos; dos(fit_Eid)];


[Ei, cur_doii] = fit_equation_linear_vbm(fit_E, fit_dos);
% Ei = linspace(min(fit_E)-3, max(fit_E)+2, 1e4);
% p = polyfit(fit_E, fit_dos, 1);
% cur_doii = polyval(p, Ei);


plot(fit_E, fit_dos, 'ok', 'linewidth', 1)
plot(Ei, cur_doii, 'r-', 'linewidth', 1)
plot(E, dos, '-', 'linewidth', 1, 'color', [0 152 255]/255)
xlim([-3 3])
[~, doi_0_id] = min(abs(cur_doii ));
loc_E_0_cbm = Ei(doi_0_id);
loc_dos_0_cbm = cur_doii(doi_0_id);

fit_E = [];
fit_dos = [];

E_pos_id = E < st_vbm & E > -3;
dos_pos = dos(E_pos_id);
E_pos = E(E_pos_id);

[~, E_pos_id] = min(abs(E_pos-st_vbm));

% [~, E_pos_id] = min(abs(dos_pos-0.001));
E_start = E_pos(E_pos_id);
fit_E_cur = linspace(E_start+vbm_step, E_start, 21);
fit_Eid = zeros(size(fit_E_cur));
for j = 1:length(fit_E_cur)
    [~, nearest_id] = min(abs(E-fit_E_cur(j)));
    fit_Eid(j) = nearest_id;
end
fit_E = [fit_E; E(fit_Eid)'];
fit_dos = [fit_dos; dos(fit_Eid)];



[Ei, cur_doii] = fit_equation_linear_vbm(-fit_E, fit_dos);
% Ei = linspace(min(fit_E)-3, max(fit_E)+2, 1e4);
% p = polyfit(fit_E, fit_dos, 1);
% cur_doii = polyval(p, Ei);

figure(1)
plot(fit_E, fit_dos, 'ok', 'linewidth', 1)
plot(-Ei, cur_doii, 'r-', 'linewidth', 1)
xlim([-5 5])
ylim([0 2000])
set(gca, 'fontname', 'arial', 'fontsize', 12, 'linewidth', 0.8, ...
    'xcolor', [0 0 0], 'ycolor', [0 0 0])

[~, doi_0_id] = min(abs(cur_doii ));
loc_E_0_vbm = -Ei(doi_0_id);
loc_dos_0_vbm = cur_doii(doi_0_id);

loc_E_0_cbm - loc_E_0_vbm

end

function [xi, yi] = fit_equation_linear_vbm(x, y)

f = fittype('A.*(x+1./E.*(5.*x.^2)+1./(E.^2).*(8.*x.^3)+1./(E.^3).*(4.*x.^4))+C', 'independent', 'x', 'coefficients', {'A', 'E', 'C'});
% f = fittype('A.*(x+(5.*x.^2)./E+(8.*x.^3)./(E.^2)+(7.*x.^3)./(E2.^2)+(4.*x.^4)./(E.^3)+(22.*x.^4)./(E*E2.^2)+(15.*x^5)./(E2.^4)+(16.*x^5)./(E^2*E2^2)+(21.*x^6)./(E*E2.^4)+(9.*x^7)./(E2.^6))+C', 'independent', 'x', 'coefficients', {'A', 'E', 'E2', 'C'});
% f = fittype('A.*(x+(5.*x.^2)./E+(15.*x.^3)./(E.^2)+(35.*x.^4)./(E.^3)+(59.*x^5)./(E.^4)+(79.*x^6)./(E.^5)+(85.*x^7)./(E.^6)+(65.*x^8)./(E.^7)+(40.*x^9)./(E.^8)+(16.*x^10)./(E.^9))+C', 'independent', 'x', 'coefficients', {'A', 'E', 'C'});
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