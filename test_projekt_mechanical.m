clear all
close all


[ndm, ndf, nnp, nel, nen, x, elem, matparam, drlt, neum, loadcurve, b] = input_Rocket_Engine_Mechanical();

plotTriaMesh(elem,x);

% separate dof types
allDofs = (1:1:nnp*ndf)';

% dofs with Dirichlet boundary conditions
numDrltDofs = ndf*length(drlt.nodes);
drltDofs = zeros(numDrltDofs,1);
for i=1:numDrltDofs/2
    node = drlt.nodes(i);
    %ldof = drlt.ldof;
    drltDofs(2*i-1) = node*ndf - 1;
    drltDofs(2*i) = node*ndf;
end

% free dofs
freeDofs = setdiff(allDofs, drltDofs);

% dofs with Neumann boundary conditions
numNeumDofs = ndf*length(neum.nodes);
neumDofs = zeros(numNeumDofs,1);
%neumNodes = [neum(1).nodes,neum(2).nodes,neum(3).nodes,neum(4).nodes,neum(5).nodes];
for i=1:numNeumDofs/2
    node = neum.nodes(i);
    %ldof = neum(i,3);
    neumDofs(2*i-1) = node*ndf - 1;
    neumDofs(2*i) = node*ndf;
end

%-----calculate normal vectors along the boundary F-I for each element----- 

%Initialise matrix containing normal coordinates for each element along FI.
%[x-coord, y-coord]
normalVectors_FI = zeros(length(neum.nodes), ndm);

%Initialise also matrix containing coordinates for the vector along the
%boundary. Who knows, we might need it.
boundaryVectors_FI = zeros(length(neum.nodes), ndm);

for i=1:length(neum.elements)
    x1=0;
    x2=0;
    y1=0;
    y2=0;
    
    % coordinates of the element nodes
    ex = x(elem(neum.nodes(i)).cn,1)*0.001;   % x-coord
    ey = x(elem(neum.nodes(i)).cn,2)*0.001;   % y-coord
    
    %wnob = which nodes on the boundary
    wnob = ismember(neum.elements(i), neum.nodes);
    
      if(wnob == [1,1,0])
         x1=ex(1);
         x2=ex(2);
         y1=ey(1);
         y2=ey(2);               
      elseif(wnob == [1,0,1])
         x1=ex(1);
         x2=ex(3);
         y1=ey(1);
         y2=ey(3);              
      elseif(wnob == [0,1,1])
         x1=ex(2);
         x2=ex(3);
         y1=ey(2);
         y2=ey(3);               
      end 
      
      %vector along the element boundary
      vec = norm([x2-x1, y2-y1]);
      
      %Fill boundaryVectors_FI with coordinates
      boundaryVectors_FI(i,1) = vec(1);
      boundaryVectors_FI(i,2) = vec(2);
      
      %Create normalvector to vector vec above.
end
%=========================================================================
% fe-analysis
%=========================================================================

% unknown node displacements
a = zeros(ndf*nnp,1);

% simulation time and time/load step
stopTime = 500.0;
dt = 10.0;
% numerical tolerance to detect equilibrium
TOL = 1e-10;

% loop over all time/load steps
time = 0;
step = 0;
while (time < stopTime)
    time = time + dt;
    step = step + 1;
    
    %Initialise global fsur (traction forces)
    fsur = zeros(ndf*nnp,1);
    
    %Initialise global fload (load forces, t.ex gravity)
    fload = zeros(ndf*nnp,1);
    
    % initialise global tangent stiffness matrix
    K = zeros(ndf*nnp,ndf*nnp);
    
    % initialise global residuum
    rsd = zeros(ndf*nnp,1);
    
    fprintf(1,'step= %2d time= %8.4e dt= %8.4e\n', step, time, dt);
    % Neumann loads for this timestep
    numNeumDofs = size(neum,1);
    neumValue = zeros(numNeumDofs,1);
    for i=1:numNeumDofs
        lcID = neum.loadcurveID;
        scale = neum.pmax * loadcurve(lcID).scalefactor;
        neumValue(i) = scale*interp1(loadcurve(lcID).time, loadcurve(lcID).value, ...
            time, 'linear', 'extrap');
    end;
    fsur(neumDofs) = neumValue;
    
    % Dirichlet boundary conditions
    drltValue = zeros(numDrltDofs,1);
    for i=1:numDrltDofs
        lcID = drlt.loadcurveID;
        scale = drlt.scale * loadcurve(lcID).scalefactor;
        drltValue(i) = scale*interp1(loadcurve(lcID).time, loadcurve(lcID).value, ...
            time, 'linear', 'extrap');
    end;
    drltValueIncr = drltValue - a(drltDofs);
 
    % element loop
    for e=1:nel
        
        % dofs belonging to element e
        edof = elem(e).edof;
        
        % coordinates of the element nodes
        ex = x((elem(e).cn),1)*0.001;   % x-coord
        ey = x((elem(e).cn),2)*0.001;   % y-coord
        
        % displacement values at the element nodes
        ae = a(edof);
        
        % area of the element
        A=1/2*det([ones(3,1) ex ey]);
              
        % B-matrix for a triangular element with 3 nodes
        B = [ ey(2)-ey(3)  0  ey(3)-ey(1)  0  ey(1)-ey(2)  0;
              0  ex(3)-ex(2)  0  ex(1)-ex(3)  0  ex(2)-ex(1);
              ex(3)-ex(2) ey(2)-ey(3) ex(1)-ex(3) ey(3)-ey(1) ex(2)-ex(1) ey(1)-ey(2) ]* 1/(2*A);
            
        % strain
        et = B*ae;
        
        % linear elasticity
        ep(1) = 1;    % plane stress
        % ep(1) = 2;    % plane strain
        [es, D, stateVar] = linElast_planeStress(et,ep,matparam, elem(e).stateVar0(1,:));
        
        % thickness of the plate
        t = matparam(4);
        
        %------Skapa och fyll fsure beroende på neum villkoren--------------
        
        %Initialise local fsur
        fsure = zeros(nen,1);
        
        %check if element is on boundary F-I
        Boolneum = ismember(elem(e).cn, neum.nodes);
        if sum(BoolRobin) == 2
        end
        
        
        % element stiffness matrix
        Ke = B'*D*B*A*t;
        
        % assemble Ke into the global stiffness matrix K
        K(edof,edof) = K(edof,edof) + Ke;
        
    end % element loop
    
    % calculate residuum and its norm
    rsd(freeDofs) =  fsur(freeDofs) - K(freeDofs,freeDofs)*a(freeDofs); 
    
    % modified right hand side
    rsd(freeDofs) = rsd(freeDofs) - K(freeDofs,drltDofs) * drltValueIncr;
    rsn = norm(rsd(freeDofs));
    fprintf(1,'residuum norm %e \n',rsn)
    if rsn > 1e-8	
    % calculate increment of node displacements
       da = zeros(size(a)); 
       da(freeDofs) = K(freeDofs,freeDofs) \ rsd(freeDofs);
       a(freeDofs) = a(freeDofs) + da(freeDofs);
       a(drltDofs) = a(drltDofs) + drltValueIncr;
       drltValueIncr = zeros(size(drltValueIncr));
    end
    % save state variables for next time/load step
    for e = 1:nel
        elem(e).stateVar0 = elem(e).stateVar;
    end;
    
    postprocessingMovie(step, time, ndm, elem, x, drltValue, drltDofs, fsur, neumDofs, a);
                      
end;