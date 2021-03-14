model=femesh('testquad4');

def.def = [0 0 1 0 0 0 0 0 2 0 0 0 0 0 3 0 0 0 0 0 4 0 0 0 ]'*[1 2]; 
def.DOF=reshape(repmat((1:4),6,1)+repmat((1:6)'/100,1,4),[],1);
def.lab={'NodeData1','NodeData2'};
sigma.EltId=[1];
sigma.data=[1 2];
sigma.lab={'EltData1' 'EltData2'};
of2vtk('fic1',model,def,sigma);

% Export to vtk with structures (.Name,.Data) 
% defined at nodes (NodeData) and at elements (EltData)
NodeData1.Name='NodeData1';NodeData1.Data=[1 ; 2 ; 3 ; 4];
NodeData2.Name='NodeData2';NodeData2.Data=[0 0 1;0 0 2;0 0 3;0 0 4];
EltData1.Name ='EltData1' ;EltData1.Data =[ 1 ];
EltData2.Name ='EltData2' ;EltData2.Data =[ 1 2 3];
of2vtk('fic2',model,NodeData1,NodeData2,EltData1,EltData2);

% Combined export
of2vtk('fic3',model,def,EltData1,EltData2);
of2vtk('fic4',model,NodeData1,NodeData2,sigma);
