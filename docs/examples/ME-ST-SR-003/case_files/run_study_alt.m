clear
close all
% Number of elements to be studied
numelemvec = [2,4,8,16];
numelems = length(numelemvec);
% Element types to be studied
% Column 1: 1 (quadrilateral), 2 (triangular)
% Column 2: 1 (linear), 2 (quadratic)
% Column 3: 0 (complete), 1 (incomplete)
elemtype = [2 1 0;
    2 2 0;
    1 1 0;
    1 2 1;
    1 2 0];

numelemtype = length(elemtype(:,1));
W=zeros(numelemtype,numelems);
NC=zeros(numelemtype,numelems,2*numelems+1);
AC=zeros(numelemtype,numelems,2*numelems+1);
WC=zeros(numelemtype,numelems,2*numelems+1);
NC=zeros(numelemtype,numelems,2*numelems+1);
AC=zeros(numelemtype,numelems,2*numelems+1);
NXC=zeros(numelemtype,numelems,2*numelems+1);
MYC=zeros(numelemtype,numelems,2*numelems+1);
QYC=zeros(numelemtype,numelems,2*numelems+1);
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
            system('gmsh -2 -order 1 cylindrical_shell.geo -o cylindrical_shell.msh -format msh2');
        else
            system('gmsh -2 -order 2 cylindrical_shell.geo -o cylindrical_shell.msh -format msh2');
        end
        system('multifebe -i cylindrical_shell_alt.dat');

        data=readtable('cylindrical_shell_alt.dat.nso', 'NumHeaderLines', 16,'FileType', 'text');
        data=table2array(data);
        W(ke,kne) = data(data(:,9)==1, 15);

        data = readtable('cylindrical_shell_alt.dat.nso', 'NumHeaderLines', 16,'FileType', 'text');
        data=table2array(data);
        data_filtered = data(data(:,11) < 0.01, [10, 15]);% x2<0.01, 24.99 x1, uz
        data_sorted = sortrows(data_filtered, 1);

        NC(ke,kne)=size(data_sorted, 1);
        AC(ke,kne,1:size(data_sorted,1))=90-acosd(data_sorted(:,1)/25);
        WC(ke,kne,1:size(data_sorted,1))=data_sorted(:,2);

        data =readtable('cylindrical_shell_alt.dat.eso', 'NumHeaderLines', 21,'FileType', 'text');
        data=table2array(data);
        data_filtered = data(data(:,16) < 0.01, [15, 9, 18, 22, 25]);% x2<0.01, x1, Element ID, Nx, My, Vy
        data_sorted = sortrows(data_filtered, 1);
        tmp = data_sorted;

        % Take only element nodes of elements with an edge on y=0
        elem=unique(tmp(:,2));
        for kel=1:length(elem)
            kevec=find(tmp(:,2)==elem(kel));
            if length(kevec)==1
                tmp(kevec,2)=0;
            end
        end
        kevec=find(tmp(:,2)~=0);
        tmp=tmp(kevec,:);
        % Calculate CDG of each element edge
        elem=unique(tmp(:,2));
        x_cdg=zeros(length(elem),1);
        for kel=1:length(elem)
            kevec=find(tmp(:,2)==elem(kel));
            x_cdg(kel)=sum(tmp(kevec,1))/length(kevec);
        end
        % Sort by element and node position
        [~,kevec]=sort(x_cdg);
        kk=1;
        tmp2=zeros(length(tmp),5);
        for kel=1:length(elem)
            knodes=find(tmp(:,2)==elem(kevec(kel)));
            [~,knsort]=sort(tmp(knodes,1));
            tmp2(kk:kk+length(knodes)-1,:)=tmp(knodes(knsort),:);
            kk=kk+length(knodes);
        end
        tmp=tmp2(1:kk-1,:);
        ENC(ke,kne)=size(tmp,1);
        EAC(ke,kne,1:size(tmp,1))=90-acosd(tmp(:,1)/25);
        NXC(ke,kne,1:size(tmp,1))=tmp(:,3);
        MYC(ke,kne,1:size(tmp,1))=tmp(:,4);
        QYC(ke,kne,1:size(tmp,1))=tmp(:,5);
    end
end
%% PLOT MAXIMUM VERTICAL DISPLACEMENT
fs=13.5;
figure
hold on
plot(numelemvec,W(1,:),'o-')
plot(numelemvec,W(2,:),'o-')
plot(numelemvec,W(3,:),'o-')
plot(numelemvec,W(4,:),'o-')
plot(numelemvec,W(5,:),'o-')
plot([2 numelemvec(end)],[-0.3024 -0.3024],'k-','LineWidth',2)
legend ('MITC3','MITC6a','MITC4','MITC8*','MITC9', 'Reference (-0.3024 in)','Location','southoutside','Orientation','Horizontal','NumColumns',3,'Interpreter','latex');
xlabel('N [-]','Interpreter','latex')
ylabel('$\mathrm{w_{midside} [in]}$','Interpreter','latex')
ylim([-0.35 -0.2])
ax = gca;
ax.FontSize = fs;
set(gcf,'Position',[100 100 1000 800])
exportgraphics (gcf,'w_midside.pdf');


%% PLOT STRESS RESULTANTS & DEFORMED SHAPE AT y=0
for ke=1:numelemtype
    figure
    if elemtype(ke,1)==2 && elemtype(ke,2)==1
        titlestr = 'MITC3';
    elseif elemtype(ke,1)==2 && elemtype(ke,2)==2
        titlestr = 'MITC6a';
    elseif elemtype(ke,1)==1 && elemtype(ke,2)==1
        titlestr = 'MITC4';
    elseif elemtype(ke,1)==1 && elemtype(ke,2)==2 && elemtype(ke,3)==1
        titlestr = 'MITC8*';
    elseif elemtype(ke,1)==1 && elemtype(ke,2)==2 && elemtype(ke,3)==0
        titlestr = 'MITC9';
    else
        error('not supported')
    end
    t=tiledlayout(1,4);
    for kne=1:numelems
        nexttile(1)
        hold on
        data=1:NC(ke,kne);
        plot(squeeze(AC(ke,kne,data)),squeeze(WC(ke,kne,data)),'o-')
        ylabel('w at y=0 in [in]','Interpreter','latex')
        ax = gca;
        ax.FontSize = fs;

        data=1:ENC(ke,kne);
        nexttile(2)
        hold on
        plot(squeeze(EAC(ke,kne,data)),squeeze(NXC(ke,kne,data)),'o-')
        ylabel('$\mathrm{N_x \ [lbf/in]}$','Interpreter', 'latex');
        ax = gca;
        ax.FontSize = fs;

        nexttile(3)
        hold on
        plot(squeeze(EAC(ke,kne,data)),squeeze(MYC(ke,kne,data)),'o-')
        ylabel('$\mathrm{M_y \ [lbf\cdot in/in]}$','Interpreter','latex')
        ax = gca;
        ax.FontSize = fs;

        nexttile(4)
        hold on
        plot(squeeze(EAC(ke,kne,data)),squeeze(QYC(ke,kne,data)),'o-')
        ylabel ('$\mathrm{Q_y \ [lbf/in]}$','Interpreter','latex')
        ax = gca;
        ax.FontSize = fs;
    end

    nexttile(1)
    plot(squeeze(AC(5,4,1:NC(5,4))),squeeze(WC(5,4,1:NC(5,4))),'k-','LineWidth',2)
    lg=legend('2 elements','4 elements','8 elements','16 elements','MITC9 16 elem.','Orientation','Horizontal','NumColumns',6,'Interpreter','latex','Fontsize',fs);
    lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout
    nexttile(2)
    plot(squeeze(EAC(5,4,1:ENC(5,4))),squeeze(NXC(5,4,1:ENC(5,4))),'k-','LineWidth',2)
    nexttile(3)
    plot(squeeze(EAC(5,4,1:ENC(5,4))),squeeze(MYC(5,4,1:ENC(5,4))),'k-','LineWidth',2)
    nexttile(4)
    plot(squeeze(EAC(5,4,1:ENC(5,4))),squeeze(QYC(5,4,1:ENC(5,4))),'k-','LineWidth',2)

    xlabel(t,'angle [$^{\circ}$]','Interpreter','latex','Fontsize',fs)
    axes('visible', 'off')
    sgtitle(titlestr,'Interpreter','latex','Fontsize',fs);

    hold off

    if strcmpi(titlestr,'MITC8*')
        titlestr='MITC8';
    end

    set(gcf,'Position',[100 100 1600 600])
    exportgraphics(gcf,strcat('stress_resultants_&_deformed_shape_',titlestr,'.pdf'));
    close all
end