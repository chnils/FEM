function [ndm, ndf, nnp, nel, nen, x, elem, matparam, drlt, neum, loadcurve, b] = input_Rocket_Engine_Mechanical()

x = dlmread('Rocket_Engine.nodes'); % 1st colum containes node number
x = x(:,2:end);              % remove 1st column

conn = dlmread('Rocket_Engine.conn'); % 1st colum containes element number
conn = conn(:,2:end);          % remove 1st column
[nel,nen] = size(conn);

% number of node points and number of spatial dimensions
[nnp, ndm] = size(x);

% number of degrees of freedom per node
ndf = ndm;  % mechanical analysis

% allocate memory for the element structure 
% 1. organize connectivity list as a struct
% 2. determine edof, element degrees of freedom
% 3. reserve memory for state variables
%    strain et [eps_xx eps_yy eps_zz gam_xy] and 
%    stress es [sig_xx sig_yy sig_zz tau_xy]
% number of state variables
nsv = 4 + 4;
elem = repmat(struct('cn',zeros(nen,1), ...
                     'edof', zeros(nen*ndf,1), ...
                     'stateVar0', zeros(1,nsv), ...
                     'stateVar', zeros(1,nsv) ...
                    ),  nel, 1 );

for e = 1:nel
    elem(e).cn = conn(e,:);
    % dofs belonging to the nodes of element e
    edof = zeros(ndf*nen,1);
    edof(1:2:ndf*nen-1) = ndf*elem(e).cn' - 1;
    edof(2:2:ndf*nen  ) = ndf*elem(e).cn';        
    elem(e).edof = edof;
end;

%----------------------------------------------------------------
% material parameter
% Young's modulus
matparam(1) = 110*10^9;
% Poisson ratio
matparam(2) = 0.36;
% density
matparam(3) = 4540;
% thickness of the plate 
matparam(4) = 0.001;

%-----------------------------------------------------------------------
% load curves
%-----------------------------------------------------------------------
id = 1;
loadcurve(id).id = id;
loadcurve(id).scalefactor =  1.0;
loadcurve(id).time =  [0.0 1.0];
loadcurve(id).value = [1.0 1.0];


%Retrieves all boundary nodes for the different boundaries.
[NodesAN,NodesAE,NodesEF,NodesFI,NodesIJ,NodesJN] = BCNodes();

% boundary conditions
% node: Node where Dirichlet condition is applied
% loadcurveID: id of load curve
% ldof: load degree of freedom (1 for heat, 1 or 2 for mechanical)
% scale: magnitude of bc

% Dirichlet boundary condition   
drlt = repmat(struct('nodes',zeros(1,length(NodesAN)),'loadcurveID',0,'ldof',0,'scale',0),1,1);
drlt.nodes = NodesAN;
drlt.loadcurveID = 1;
drlt.ldof = 0;
drlt.scale = 0;
%antar att scale (och ldof?) ska vara noll eftersom v�ggen �r or�rlig.  

%neum elements, going from F to I
neumElements = [71,242,176,168,160,144,137,136,126,112,107,124,64,193,62,99,101,152,15,...
 102,55,114,130,138,148,162,170,178,186,50,52];
neumElements = sort(neumElements);

% Neumann boundary condition  
neum = repmat(struct('nodes',zeros(1,length(NodesFI)),'pmax',0,'loadcurveID',id,'elements',neumElements),1,1); 

%F-I boundary
neum.nodes = NodesFI;
neum.pmax = 10*10^6; 


% % Neumann boundary condition   
% neum = repmat(struct('nodes',zeros(1,length(NodesFI)),'pmax',0,'loadcurve',id),5,1); 
% 
% %A-E boundary
% neum(1).nodes = NodesAE;
% neum(1).pmax = 0;
% 
% %E-F boundary
% neum(2).nodes = NodesEF;
% neum(2).pmax = 0;
% 
% %F-I boundary
% neum(3).nodes = NodesFI;
% neum(3).pmax = 10*10^6; 
% 
% %I-J boundary
% neum(4).nodes = NodesIJ;
% neum(4).pmax = 0;
% 
% %J-N boundary
% neum(5).nodes = NodesJN;
% neum(5).pmax = 0; 


% body force aka gravity, Ska g=9.81 in p� plats tv� f�r att beskriva att
% dess verkan sker i y-led.
b = [0; 9.81];