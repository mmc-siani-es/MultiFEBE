%% OCTAVE SCRIPT 

clear all
close all

% Number of elements to be studied
numelemvec = [2 4 8 16];
numelems = length(numelemvec);

% Element types to be studied
% Column 1: 1 (quadrilateral), 2 (triangular)
% Column 2: 1 (linear), 2 (quadratic)
% Column 3: 0 (complete), 1 (incomplete)
elemtype = [2 1 0;  % mitc3
            2 2 0;  % mitc6a
            1 1 0;  % mitc4
            1 2 1;  % mitc8*
            1 2 0]; % mitc9       
            
numelemtype = length(elemtype(:,1));

W=zeros(numelemtype,numelems);

for ke=1:numelemtype

  disp('Element')
  disp(ke)

  for kne=1:numelems
  
    numelem = numelemvec(kne);
    
    disp('Element')
    disp(numelem)    
    
    f=fopen('ne.geo','w');
    fprintf(f,'ne=%g;\n',numelem);
    fprintf(f,'elem=%g;\n',elemtype(ke,1));
    fprintf(f,'incom=%d;\n',elemtype(ke,3));
    fclose(f);      
    if elemtype(ke,2)==1
      system('gmsh -2 -order 1 spherical_dome.geo -o  spherical_dome.msh -format msh2');
    else
      system('gmsh -2 -order 2 spherical_dome.geo -o  spherical_dome.msh -format msh2');
    end

    system('multifebe -i spherical_dome.dat');

    data = importdata('spherical_dome.dat.nso', ' ', 16);
    W(ke,kne)=data.data(1,13);
    
  end
end

%% PLOT MAXIMUM VERTICAL DISPLACEMENT

figure
hold on
plot(numelemvec,W(1,:),'o-')
plot(numelemvec,W(2,:),'o-')
plot(numelemvec,W(3,:),'o-')
plot(numelemvec,W(4,:),'o-')
plot(numelemvec,W(5,:),'o-')
plot([2 numelemvec(end)],[0.0940 0.0940],'k-','LineWidth',2)
legend('MITC3','MITC6a','MITC4','MITC8*','MITC9','Reference (0.0924 in)','Interpreter','latex');
legend('location','south')
xlabel('N [-]','Interpreter','latex')
ylabel('$\mathrm{u_x \ in \ (10,0,0)}$','Interpreter','latex')
ylim([0.05 0.1])
ax = gca;
ax.FontSize = 13.5;
set(gcf,'Position',[100 100 1000 800])
exportgraphics (gca,'ux_at_load.pdf');

close all
