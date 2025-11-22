clear
close all
% Number of elements to be studied
numelemvec = [2 4 8 16];
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
W=zeros(numelemtype,numelems,2);
forc_app=["Concentrated","Uniform"];


new_conditions_nodes_dist = {
    '[conditions over nodes]'
    '' % Blank line for separation
    'group 1: 0 0'
    '         0 0'
    '         0 0'
    '         0 0'
    '         0 0'
    '         0 0'
    '         ' % Blank line
    'group 2: 0 0'
    '         0 0'
    '         0 0'
    '         0 0'
    '         0 0'
    '         0 0'
    '         ' % Blank line
    };

for fa=1:2 %force application
    if fa==1
        new_conditions_nodes=new_conditions_nodes_dist;
        new_conditions_nodes(17:23,:) = {
            'node 1:  0 0.'
            '         0 0.'
            '         1 1.d-4'
            '         0 0.'
            '         0 0.'
            '         0 0.'
            '         ' % Blank line
            };

        % New content for the 'conditions over fe elements' block
        new_conditions_fe = {};
    else
        new_conditions_nodes = new_conditions_nodes_dist;
        new_conditions_fe = {
            '[conditions over fe elements]'
            'group 3: global'
            '         0.d0 0.d0 1.d-4'
            };
    end

    for ke=1:numelemtype
        disp('Element')
        disp(ke)
        for kne=1:numelems
            numelem = numelemvec(kne);
            disp('Element')
            disp(numelem)
            disp(forc_app(fa))

            f=fopen('ne.geo','w');
            fprintf(f,'ne=%g;\n',numelem);
            fprintf(f,'elem=%g;\n',elemtype(ke,1));
            fprintf(f,'incom=%d;\n',elemtype(ke,3));
            fclose(f);
            modify_rectangular_plate_data(new_conditions_nodes,new_conditions_fe)

            if elemtype(ke,2)==1
                system('gmsh -2 -order 1 Rectangular_plate.geo -o Rectangular_plate.msh -format msh2');
            else
                system('gmsh -2 -order 2 Rectangular_plate.geo -o Rectangular_plate.msh -format msh2');
            end
            system('multifebe -i Rectangular_plate.dat');

            data = importdata('Rectangular_plate.dat.nso', ' ', 16);

            W(ke,kne,fa) = data.data(1,15);
            clearvars data
        end
    end
end

%% PLOT MAXIMUM VERTICAL DISPLACEMENT
for fa=1:2 %force application
    figure(fa)
    hold on
    plot(numelemvec,W(1,:,fa),'o-')
    plot(numelemvec,W(2,:,fa),'o-')
    plot(numelemvec,W(3,:,fa),'o-')
    plot(numelemvec,W(4,:,fa),'o-')
    plot(numelemvec,W(5,:,fa),'o-')
    if fa==1
        plot([2 numelemvec(end)],[5.60 5.60],'k-','LineWidth',2)
        legend ('MITC3','MITC6a','MITC4','MITC8*','MITC9', 'Reference (5.60 in)','Location','best','Interpreter','latex');
        ylim([0 6])
    else
        plot([2 numelemvec(end)],[1.26 1.26],'k-','LineWidth',2)
        legend ('MITC3','MITC6a','MITC4','MITC8*','MITC9', 'Reference (1.26 in)','Location','best','Interpreter','latex');
        ylim([0 2])
    end
    xlabel('N [-]','Interpreter','latex')
    ylabel('$\mathrm{w_{midside} [in]}$','Interpreter','latex')
    exportgraphics (gca,'w_midside_force_application'+forc_app(fa)+'.pdf');
    hold off

    close all

end

function modify_rectangular_plate_data(new_conditions_nodes, new_conditions_fe)
% File name to edit
filename = 'Rectangular_plate.dat';

% -------------------------------------------------------------------------
% 1. Read the file and find the indices of the blocks
% -------------------------------------------------------------------------
fid = fopen(filename, 'r');
if fid == -1
    error('Could not open the file.');
end
file_content = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', '');
fclose(fid);
file_content = file_content{1};

% Find the starting lines of each block
start_conditions_nodes = find(strcmp(file_content, '[conditions over nodes]'), 1);
start_conditions_fe = find(strcmp(file_content, '[conditions over fe elements]'), 1);

% -------------------------------------------------------------------------
% 2. Assemble the new file
% -------------------------------------------------------------------------

% Find the end of each block
end_of_conditions_nodes = 0;
if ~isempty(start_conditions_nodes)
    idx = start_conditions_nodes + 1;
    while idx <= length(file_content) && ~startsWith(file_content{idx}, '[')
        idx = idx + 1;
    end
    end_of_conditions_nodes = idx - 1;
end

end_of_conditions_fe = 0;
if ~isempty(start_conditions_fe)
    idx = start_conditions_fe + 1;
    while idx <= length(file_content) && ~startsWith(file_content{idx}, '[')
        idx = idx + 1;
    end
    end_of_conditions_fe = idx - 1;
end

% Build the new content step by step
new_content = {};
new_content = [new_content; file_content(1:start_conditions_nodes - 1)];

% Add the new blocks, replacing the old ones
new_content = [new_content; new_conditions_nodes];
if ~isempty(new_conditions_fe)
    new_content = [new_content; new_conditions_fe];
end

% Add the rest of the file content after the last modified block
last_block_end = end_of_conditions_nodes;

if end_of_conditions_fe > last_block_end
    last_block_end = end_of_conditions_fe;
end

if last_block_end < length(file_content)
    new_content = [new_content; file_content(last_block_end + 1:end)];
end

% -------------------------------------------------------------------------
% 3. Write the new content to the original file
% -------------------------------------------------------------------------
fid = fopen(filename, 'w'); % 'w' para sobrescribir
if fid == -1
    error('Could not write to the file.');
end
fprintf(fid, '%s\n', new_content{:});
fclose(fid);
end