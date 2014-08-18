% Ploteo de la pieza.
% Sirve para la previzualización y para la superposición final con la parte deformada.


elems =  load('ELIST.txt');   % Esto se necesita para saber los nodos que conforman un elemento.
%load('nodes.txt');            % Posicion inicial de los nodos.
 elems(:,2:6) = [];

%load('V5_2.txt'); 


    nodes_def = nodes;%zeros(size(nodes,1),4);
if pre == 0    
    nodes_def(:,1) = nodes(:,1);  %No se si este paso es necesario.
    
    scaleFact = 0.1;
    nodes_def(:,2:4) =  NodeSort(:,2:4);
end

%nodes_def(:,2:4) = nodes(:,2:4) + u_conv*scaleFact; 


%nodes_def(:,2:4) = nodes(:,2:4) + V5_2(:,3:5)*scaleFact;   % Para ser usado con los modos extraidos de Ansys.  
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
        
        node_ini = nodes(elems(i,j+5),2:4); elems(i,j+1);
        node_fin = nodes(elems(i,j+6),2:4); elems(i,j+2);
        
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
   
     %if pre == 0
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
     %end
    
        
        
end
