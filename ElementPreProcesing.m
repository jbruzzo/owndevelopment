% Preprocesing the fem element: Drawing the part and extracting the
% constraint modes.


% 
clear all
close all
clc

% Location of the nodes

% Interface node p.

% ^
% | y
%     -> z

% 2094 2095   2   4   3
  
% 1452 1453   1   6   5

% 1450 1451 968 969 967



% Interface node q.

%       ^
%       | y
% <- Z



%   9  10   8 2097 2096

%  11  12   7 1457 1456

% 970 972 971 1455 1454





 load('nodes.txt');
 Kfem = load('stiffness_K.dat');
 Mfem = load('Mass_M.dat');
 
 
 %[V,D] = eig(Kfem,Mfem);
 
 

% Detecting the nodes on the interfaces...


elem_length = max(nodes(:,2));    % Length of the element.


npcount = 1;    % Interface node counters.
nqcount = 1;

for i = 1:length(nodes(:,1))
    if nodes(i,2) == 0
        inter_p(npcount,1) = nodes(i,1);  % Vector with the p-interface nodes.
        npcount = npcount+1;
    
    elseif nodes(i,2) == elem_length
        inter_q(nqcount,1) = nodes(i,1);  % Vector with the q-interface nodes.
        nqcount = nqcount+1;
           
    end
end

%% Rearranging the stiffness matrix to comply with
%
%     Kfem cc |  Kfem ci      Uc        Fc
%     ------------------  x        =  
%     Kfem ic |  Kfem ii      Ui        0




clear Knew Kfem2

Kfem2 = Kfem;

vcount = 1;    % Interface node counters.
kcount = 1;

Knew = zeros(3*length(nodes));

for i = 1:length(nodes(:,1))
    for j = 1:length(inter_p(:,1))
        if nodes(i,1) == inter_p(j,1)
           % Knew((3*kcount-2):3*kcount,:) = Kfem((3*inter_p(j,1)-2):3*inter_p(j,1),:);
            %Knew(:,(3*kcount-2):3*kcount) = Kfem(:,(3*inter_p(j,1)-2):3*inter_p(j,1));
            VecOrder(kcount,1) = nodes(i,1);
            
            kcount = kcount + 1;
            vcount = vcount + 1;
        end
    end
end

for i = 1:length(nodes(:,1))
    for j = 1:length(inter_q(:,1))
        if nodes(i,1) == inter_q(j,1)
            %Knew((3*kcount-2):3*kcount,:) = Kfem((3*inter_q(j,1)-2):3*inter_q(j,1),:);
            %Knew(:,(3*kcount-2):3*kcount) = Kfem(:,(3*inter_q(j,1)-2):3*inter_q(j,1));
            VecOrder(kcount,1) = nodes(i,1);
            
            kcount = kcount + 1;
            vcount = vcount + 1;
        end
    end
end





NodeNum = nodes(:,1);


for i=1:length(VecOrder(:,1))
    NodeNum(VecOrder(i,1),:) = 1;
end

rcount = 1;

for i = 1:length(NodeNum(:,1))
    if  NodeNum(i,1) ~= 1
        NodeNum2(rcount,1) = NodeNum(i,1);
        rcount = rcount + 1;
    end
    
end



VecOrderFinal = [VecOrder;NodeNum2];  % Final order of the nodes arranged in the stiffness matrix.

for i=1:length(VecOrderFinal(:,1))

    nodes2(i,1:3) = nodes(VecOrderFinal(i,1),2:4);  % Final location of the nodes according to the new arrangement. 

end




vcount = 1;    % Interface node counters.
kcount = 1;

Knew1 = zeros(3*length(nodes));
Knew2 = zeros(3*length(nodes));
Knew3 = zeros(3*length(nodes));
Knew4 = zeros(3*length(nodes));




for i = 1:length(VecOrderFinal(:,1))
    
    % Reordenando todas las filas
    Knew1(3*i-2:3*i,:) = Kfem(3*VecOrderFinal(i,1)-2:3*VecOrderFinal(i,1),:);
    
end


for i = 1:length(VecOrderFinal(:,1))
    
    % Reordenando todas las columnas.
    Knew2(:,3*i-2:3*i) = Knew1(:,3*VecOrderFinal(i,1)-2:3*VecOrderFinal(i,1),:);
    
end






% This part was the one previously working.

% for i = 1:length(nodes(:,1))
%     for j = 1:length(VecOrderFinal(:,1))
%         if nodes(i,1) == VecOrderFinal(j,1)
%            Knew((3*kcount-2):3*kcount,:) = Kfem((3*VecOrderFinal(j,1)-2):3*VecOrderFinal(j,1),:);
%            Knew(:,(3*kcount-2):3*kcount) = Kfem(:,(3*VecOrderFinal(j,1)-2):3*VecOrderFinal(j,1));
%            
%             
%             kcount = kcount + 1;
%             vcount = vcount + 1;
%         end
%     end
% end


% for i = 1:length(VecOrder(:,1))
%     Kfem2(3*VecOrder(i,1)-2:3*VecOrder(i,1),:) = [];
% end
% 
K = Knew2;  % Rearranged stifness matrix.


% Constraint mode V1

% Resizing the stiffness matrix to comply with the constraints imposed.

% Eliminating the rows and columns from the rearranged stiffness matrix
% related to the constrained dof's.



K(3*length(inter_p)+1:3*(length(inter_p)+length(inter_q)),:) = []; % Rows related to the interface q.

K(:,3*length(inter_p)+1:3*(length(inter_p)+length(inter_q))) = []; % Columns related to the interface q;


% Now the fixed dof of the interface point p, i.e., y and z directions.

for i=length(inter_p):-1:1
    i;
    K(3*i-1:3*i,:) = []; % Row elimination
    K(:,3*i-1:3*i) = []; % Column elimination
    
       
    
end






ci = 3*size(VecOrder,1) - 3*length(inter_q) - 2*length(inter_p);
ii = ci+1;

Kii = K(ii:end,ii:end);
Kci = K(ci+1:end,1:ci);




rcount = 1;

uc = zeros(ci,1);


for i=1:ci
        uc(i,1) = 1;
end





ui = -Kii\(Kci*uc);


%%% Proximo paso es hacer el delta de la posicion de los nodos







% forming the full nodal pos vector

NodeFull = zeros(length(VecOrderFinal),4);

NodeFull(:,1) = VecOrderFinal;

NodeFull(:,2:4) = nodes2;


ScaleFactor = 0.07;

uc_s = uc*ScaleFactor;
ui_s = ui*ScaleFactor;


% Expanding the local displacements

uc_exp = zeros(2*3*length(inter_p),1);  % considering two sets of interface nodes.

for i = 1:length(inter_p)
    
    uc_exp(3*i-2,1) = uc_s(i,1);
        
end

u = [uc_exp;ui_s];

ucount = 1;      % Converting the list u into coordinate wise form.
for i=1:size(nodes,1)
    for j = 1:3
        u_conv(i,j) = u(ucount,1);
        
        ucount = ucount +1;
    end
    
end

NodeFull(:,2:4) = NodeFull(:,2:4) + u_conv;



    



% sorting the elements of matrix containing the new position of the node.
% It has to be in the same format of "nodes", so it can be printed.


for i = 1:length(VecOrderFinal)
    for j = 1:length(VecOrderFinal)
        if i == VecOrderFinal(j,1)
            NodeSort(i,:) = NodeFull(j,:);
            
        end
        
    end
    
end



%%

% % % clear all
% % % close all
% % % clc
% % % 
% % % 
elems =  load('ELIST.txt');
load('nodes.txt'); 
 elems(:,2:6) = [];

%load('V2_2.txt'); 

nodes_def = zeros(size(nodes,1),4);

nodes_def(:,1) = nodes(:,1);



scaleFact = 0.1;


nodes_def(:,2:4) = NodeSort(:,2:4); 
%nodes_def(:,2:4) = nodes(:,2:4) + V2_2(:,3:5)*scaleFact;   % Para ser usado con los modos extraidos de Ansys.  
%nodes = NodeSort; 
% %%
figure(1)
% 
for i = 1:size(elems,1)
    for j=1:3
        %plotting the faces
        node_ini = nodes(elems(i,j+1),2:4);
        node_fin = nodes(elems(i,j+2),2:4);
        
        plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)])
        hold on
        
        node_ini = nodes(elems(i,j+5),2:4); elems(i,j+1)
        node_fin = nodes(elems(i,j+6),2:4); elems(i,j+2)
        
        plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)])
        
        %plotting the longitudinal lines
        node_ini = nodes(elems(i,j+1),2:4);% elems(i,j+1)
        node_fin = nodes(elems(i,j+5),2:4);% elems(i,j+2)
        
        plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)])
  
        axis([-0 2 -0.5 0.5 -0.5 0.5])
        
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        

    end
     
    for j=1:3
        %plotting the faces
        node_ini = nodes_def(elems(i,j+1),2:4);
        node_fin = nodes_def(elems(i,j+2),2:4);
        
        plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)],'r')
        hold on
        
        node_ini = nodes_def(elems(i,j+5),2:4); %elems(i,j+1)
        node_fin = nodes_def(elems(i,j+6),2:4);% elems(i,j+2)
        
        plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)],'r')
        
        %plotting the longitudinal lines
        node_ini = nodes_def(elems(i,j+1),2:4);% elems(i,j+1)
        node_fin = nodes_def(elems(i,j+5),2:4);% elems(i,j+2)
        
        plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)],'r')
  
        axis([-0 2 -0.5 0.5 -0.5 0.5])
        

     end
    
        
        
end
%%


% % % %% Mapping back from the boundary nodes.
% % % 
% % % clear all
% % % clc
% % % 
% % % 
% % % y0 = [  0, 0, 0,  1,  0,  0,  0, 2,  0,  0,   1,   0,   0,   0];
% % % 
% % % 
% % % e0p = y0(4);
% % % e1p = y0(5);
% % % e2p = y0(6);
% % % e3p = y0(7);
% % % 
% % % e0q = y0(11);
% % % e1q = y0(12);
% % % e2q = y0(13);
% % % e3q = y0(14);
% % % 
% % % 
% % % % Boundary node p
% % % 
% % % ep = [e1p;e2p;e3p];
% % % 
% % % Deltap = [-ep,e0p*eye(3) + skew(ep)];     % Eq. 14
% % % 
% % % % Boundary node q
% % % 
% % % eq = [e1q;e2q;e3q];
% % % 
% % % Deltaq = [-eq,e0q*eye(3) + skew(eq)];
% % % 
% % % 
% % % TransM(1:3,1:3) = (eye(3));
% % % TransM(4:6,4:7) = 2*Deltap;%Re.'*Rp.'*Deltap;
% % % TransM(7:9,8:10) = (eye(3));
% % % TransM(10:12,11:14) = 2*Deltaq;%Re.'*Rq.'*Deltaq;
% % % 
% % % 
% % % 
% % %     % Set with just 3 DOF per node.
% % %     
% % %     load('V1_2.txt');
% % %     load('V2_2.txt');
% % %     load('V3_2.txt');
% % %     load('V4_2.txt');
% % %     load('V5_2.txt');
% % %     load('V6_2.txt');
% % %     load('V7_2.txt');
% % %     load('V8_2.txt');
% % %     load('V9_2.txt');
% % %     load('V10_2.txt');
% % %     load('V11_2.txt');
% % %     load('V12_2.txt');
% % %     
% % % 
% % %      for i=1:size(V1_2,1)       
% % %         V1arr((3*i-2):3*i,1)  = [V1_2(i,3);V1_2(i,4);V1_2(i,5)];
% % %         V2arr((3*i-2):3*i,1)  = [V2_2(i,3);V2_2(i,4);V2_2(i,5)];
% % %         V3arr((3*i-2):3*i,1)  = [V3_2(i,3);V3_2(i,4);V3_2(i,5)];
% % %         V4arr((3*i-2):3*i,1)  = [V4_2(i,3);V4_2(i,4);V4_2(i,5)];
% % %         V8arr((3*i-2):3*i,1)  = [V8_2(i,3);V8_2(i,4);V8_2(i,5)];
% % %         V9arr((3*i-2):3*i,1)  = [V9_2(i,3);V9_2(i,4);V9_2(i,5)];
% % %         
% % %         V5arr((3*i-2):3*i,1)  = [V5_2(i,3);V5_2(i,4);V5_2(i,5)];
% % %         V6arr((3*i-2):3*i,1)  = [V6_2(i,3);V6_2(i,4);V6_2(i,5)];
% % %         V7arr((3*i-2):3*i,1)  = [V7_2(i,3);V7_2(i,4);V7_2(i,5)];
% % %         V10arr((3*i-2):3*i,1) = [V10_2(i,3);V10_2(i,4);V10_2(i,5)];
% % %         V11arr((3*i-2):3*i,1) = [V11_2(i,3);V11_2(i,4);V11_2(i,5)];
% % %         V12arr((3*i-2):3*i,1) = [V12_2(i,3);V12_2(i,4);V12_2(i,5)];
% % %         
% % %      end
% % %      
% % % % Matrix of eigenvectors for the mass matrix
% % % 
% % % V = [V1arr,V2arr,V3arr,V4arr,V5arr,V6arr,V7arr,V8arr,V9arr,V10arr,V11arr,V12arr];  % from Eq. 2
% % % 
% % % 
% % % 
% % % u = V*(TransM*y0');
% % % 
% % % 
% % % 
% % % ucount = 1;      % Converting the list u into coordinate wise form.
% % % for i=1:size(V,1)/3
% % %     for j = 1:3
% % %         u_conv(i,j) = u(ucount,1);
% % %         
% % %         ucount = ucount +1;
% % %     end
% % %     
% % % end
% % % 
% % % 
% % % elems =  load('ELIST.txt');
% % % 
% % %  elems(:,2:6) = [];
% % % nodes00 = load('nodes.txt'); 
% % %  
% % % 
% % % nodes(:,1) = nodes00(:,1);
% % % 
% % % nodes(:,2:4) = u_conv(:,1:3);
% % % 
% % % scaleFact = 0.1;
% % % 
% % % nodes_def(:,2:4) = nodes00(:,2:4) + u_conv*scaleFact;  
% % %  
% % % figure(1)
% % % % 
% % % for i = 1:size(elems,1)
% % %     for j=1:3
% % %         %plotting the faces
% % %         node_ini = nodes_def(elems(i,j+1),2:4);
% % %         node_fin = nodes_def(elems(i,j+2),2:4);
% % %         
% % %         plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)])
% % %         hold on
% % %         
% % %         node_ini = nodes_def(elems(i,j+5),2:4); elems(i,j+1)
% % %         node_fin = nodes_def(elems(i,j+6),2:4); elems(i,j+2)
% % %         
% % %         plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)])
% % %         
% % %         %plotting the longitudinal lines
% % %         node_ini = nodes_def(elems(i,j+1),2:4);% elems(i,j+1)
% % %         node_fin = nodes_def(elems(i,j+5),2:4);% elems(i,j+2)
% % %         
% % %         plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)])
% % %   
% % %         axis([-0 2 -0.5 0.5 -0.5 0.5])
% % %         
% % %         xlabel('X')
% % %         ylabel('Y')
% % %         zlabel('Z')
% % %         
% % % 
% % %     end
% % %      
% % % %     for j=1:3
% % % %         %plotting the faces
% % % %         node_ini = nodes_def(elems(i,j+1),2:4);
% % % %         node_fin = nodes_def(elems(i,j+2),2:4);
% % % %         
% % % %         plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)],'r')
% % % %         hold on
% % % %         
% % % %         node_ini = nodes_def(elems(i,j+5),2:4); %elems(i,j+1)
% % % %         node_fin = nodes_def(elems(i,j+6),2:4);% elems(i,j+2)
% % % %         
% % % %         plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)],'r')
% % % %         
% % % %         %plotting the longitudinal lines
% % % %         node_ini = nodes_def(elems(i,j+1),2:4);% elems(i,j+1)
% % % %         node_fin = nodes_def(elems(i,j+5),2:4);% elems(i,j+2)
% % % %         
% % % %         plot3([node_ini(1,1) node_fin(1,1)],[node_ini(1,2) node_fin(1,2)],[node_ini(1,3) node_fin(1,3)],'r')
% % % %   
% % % %         axis([-0 2 -0.5 0.5 -0.5 0.5])
% % % %         
% % % % 
% % % %      end
% % %     
% % %         
% % %         
% % % end
