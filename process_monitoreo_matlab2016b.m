%SCRIPT PARA MOSTRAR LA EVOLUCION DE LA TMIN
%TOMANDO COMO REFERENCA LA CLIMATOLOGIA EN HUAYAO
%MODIFICAR: DIR_data: directorio del archivo de datos Tmin
%           namefile: nombre del archivo de datos con columnas
%			anho | mes | dia | Tmin

close all
clear all
clc 

set(0,'defaulttextfontsize',10)
set(0,'defaultaxesfontsize',9)

%******************************
%**	Directorio de datos
DIR_data = './';
%******************************
%**	Nombre del archivo de datos
namefile = 'Tn_Huayao_90_19.txt';

%************************************
%***	Cargando data
%************************************
data = load([DIR_data '/' namefile]);
t_vec = data(:,1:3);
t_num = datenum(t_vec);
 Tmin = data(:,4);

%************************************
%***	Creando anho climatologico
%************************************
t_num_clim = datenum([2014 01 01 0 0 0]):datenum([2014 12 31 0 0 0]);
t_vec_clim = datevec(t_num_clim);
    i_tick = find(t_vec_clim(:,3)==01);
     xtick = t_num_clim(i_tick);
xticklabel = datestr(xtick,'mm');

%************************************
%***	buscando todos los anhos
yyyys = unique(t_vec(:,1));
Ny = length(yyyys);
%************************************
%***	Clim para c/dia
for i = 1:365
    t_vec_ref = t_vec_clim(i,:);
    X=[];
  for yyyy = yyyys'
    %*** Specific day to evaluate
    t_vec_ref = [yyyy t_vec_ref(2:6)];
    t_num_ref = datenum(t_vec_ref);
    %*** Window around the specific day
    i_yy_day = find(t_num>=t_num_ref-15 & t_num<=t_num_ref+15); %-15->+15 dias
    %*** Limit of years for climatology
    i_yy_lim = find(t_num>=datenum([1990 01 01 00 00 00]) & t_num<datenum([2020 01 01 00 00 00]));
    i_yy     = intersect(i_yy_day,i_yy_lim);
    X = [X; Tmin(i_yy)]; 
  end
    Pr(1:5,i)  = prctile(X,[05 25 50 75 95]);
end
%***	Guardando percentiles
	fid = fopen('Tmin_percentiles_Tmin.txt','w')
	t_clim = t_vec_clim; t_clim(:,1) = 2099;
	fprintf(fid,'%4.4i %2.2i %2.2i  %6.1f %6.1f %6.1f %6.1f %6.1f \n',[t_clim(:,1:3) Pr']')
	fclose(fid)
clear t_clim

%***************************************************
%***	Completando dias q faltan en el year con NaN
%***************************************************
t_vec_x = t_vec(end,:);
t_vec_x(2:6) = [12 31 00 00 00] ;
t_num_x = datenum(t_vec_x);
t_num_N = t_num(1):t_num_x;
t_vec_N = datevec(t_num_N);
%***	Detectando y eliminando 29Feb's
i_x = find(t_vec_N(:,2) == 02 & t_vec_N(:,3) == 29);
t_num_N(i_x) = [];
t_vec_N = datevec(t_num_N);
%***	Associating data
Tmin_N = t_num_N*NaN;
[x i_x ii] = intersect(t_num_N,t_num);
Tmin_N(i_x) = Tmin(ii);

%************************************%************************************
%***	FIGURAS CLIM
%************************************%************************************
fig_clim = figure('Visible','Off')
DX = 16;
DY = 8;
papersize = [DX DY];
  set(fig_clim,'paperunits','Centimeters')
  set(fig_clim,'paperorientation','landscape')
  set(fig_clim,'papersize',papersize)
  set(fig_clim,'paperposition',[0 0 papersize])
  pos = [1.2/DX 1.2/DY 14/DX 6/DY]
subplot('position',pos)
hold on
h05_95 = patch([t_num_clim t_num_clim(end:-1:1)],[Pr(1,:) Pr(5,end:-1:1)],'c')
h25_75 = patch([t_num_clim t_num_clim(end:-1:1)],[Pr(2,:) Pr(4,end:-1:1)],'g')
   h50 = plot([t_num_clim],[Pr(3,:)],'k','linewidth',2)
  set(h05_95,'facecolor',[1 1 1]*0.90,'edgecolor',[1 1 1]*0.75)
  set(h25_75,'facecolor',[1 1 1]*0.75,'edgecolor',[1 1 1]*0.50)
  box on
  ylim([-6 11])
  xlim([t_num_clim(1)-001 t_num_clim(end)+1])
  set(gca,'ytick',[-6:3:12])
  set(gca,'xtick',xtick)
  set(gca,'xticklabel',cellstr(xticklabel))
  ylabel('Tmin ({\circ}C)')
  xlabel('Meses')
hzero = plot(t_num_clim,Pr(3,:)*0,'-b','linewidth',1)
hnow = plot(t_num_clim,Tmin_N(end-364:end),'.-r','linewidth',1)
  leg=legend([h05_95 h25_75 h50 hnow],'Perc. [5-95]','Perc. [25-75]','Mediana','2021')
  legend boxoff
  set(leg,'location','southwest')

print('Fig_Tmin_climate_monitoreo','-dpng','-r500')
%print('Fig_Tmin_climate_monitoreo','-dpdf')				%GENERA UNA FIGURA EN pdf

%**********************************************
%*** Calculando Promedios mensuales	
%*** Climatologia para el promedio
i_y = 0; %cont de years
i_m = 0; %cont de meses
Tmin_m = [];
Tmin_m_clim = [];
for y=t_vec_N(1,1):t_vec_N(end,1)
  i_y = i_y+1;
  for m=1:12
	ind_x = find(t_vec_N(:,1)==y & t_vec_N(:,2)==m);
	i_m = i_m+1;
	Tmin_m(i_m,1) = nanmean(Tmin_N(ind_x));
	if(y<=2019)
	  Tmin_m_clim(i_y,m) = nanmean(Tmin_N(ind_x));
	end
	t_vec_m(i_m,:) = t_vec_N(ind_x(1),1:2);
  end
end
Tmin_m_clim = mean(Tmin_m_clim)';
Tmin_m_clim = repmat(Tmin_m_clim,i_y,1);
%***    Guardando climatologías de promedios mensuales
    fid = fopen('Tmin_clim_prom_mensual.txt','w')
    t_clim = t_vec_clim; t_clim(:,1) = 2099;
    fprintf(fid,'Year  mm  clim  avg anom \n')
    fprintf(fid,'%4.4i %2.2i  %6.1f %6.1f %6.1f \n',[t_vec_m Tmin_m_clim Tmin_m Tmin_m-Tmin_m_clim]')
    fclose(fid)
return

