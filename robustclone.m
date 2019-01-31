tic
D=csvread('example.csv',1,1);
[A1,E1]=exact_alm_rpca(D);  % A1 is the required low rank matrix


AA1=int8(A1);  % AA1 is the integerized A1
EE1=int8(E1);
%save('AA1.mat','AA1')


length(find(AA1==2))
[AA1_1,order1]=unique(AA1,'rows');
AA1_1(all(AA1_1==0,2),:) = [];
[AA1_2, order2]= unique(AA1_1','rows');
AA1_2=AA1_2';
[p,q]=size(AA1_2);
AA1_3=AA1_2(:,q:-1:1);  %AA1_3 is the simplified and sorted matrix of AA1
figure(1)
heatmap(AA1_3)


clone_genes={}; % the mutation indexes contained in each mutation class 
clone_cells={};  %  the cell indexes contained in each clone
[nc1,ng1]=size(AA1_1);
[nc2,ng2]=size(AA1_3);
for i=1:ng2
    beclone=zeros(1,ng1);
    ge=AA1_3(:,i);
    for j=1:ng1
        cge=AA1_1(:,j);
        beclone(1,j)=pdist2(double(ge'),double(cge'),'hamming');
    end
    clone_ge=find(beclone==0);
    clone_genes{i}= clone_ge;
end


clone_mutation={}; %the matation sites for each clone
for i=1:nc2
    ce=AA1_3(i,:);
    ce_site=find(ce~=0);
    if length(ce_site)>0
        clone_mutation{i}= clone_genes{ce_site(1)};
        if length(ce_site)>1
            for j=2:length(ce_site)
                clone_mutation{i}= [clone_mutation{i}, clone_genes{ce_site(j)}];
            end
        end
    end
end


[nc3,ng3]=size(AA1);
for i=1:nc1
    beclone=zeros(1,nc3);
    ce=AA1(order1(i),:);
    for j=1:nc3
        cce=AA1(j,:);
        beclone(1,j)=pdist2(double(ce),double(cce),'hamming');
    end
    clone_ce=find(beclone==0);
    clone_cells{i}= clone_ce;
end

 % Huffman tree
k=0;
code={};
leaves={};
for i=1:p
    cell=AA1_2(i,:);
    t = cell(find(cell, 1, 'first'):end);
    tt = t(length(t):-1:1);
    leaves{i}=tt;
    sorted_str={};
    for j = 1:length(tt)
        sorted_str{j} = strrep(num2str(tt(1:j)), ' ', '');
    end
    k=k+length(sorted_str);
    code{i} = sorted_str;
end

IDS={};
for i=1:length(code)
    for j=1:length(code{i})
        IDS=[IDS,code{i}(j)];
    end
end

uIDS=unique(IDS);
adjmatrix=zeros(length(uIDS),length(uIDS));
for i=1:length(code)
    for j=1:(length(code{i})-1)
        scode1=code{i}{j};
        scode2=code{i}{j+1};
        [bool,p1]=ismember(scode1,uIDS);
        [bool,p2]=ismember(scode2,uIDS);
        adjmatrix(p1,p2)=1;
    end
end

label={};
label{1}='root';
for i=2:length(uIDS)
     label{i}=strrep(['N',num2str(i)], ' ', '');
end

clone_position=zeros(1,length(leaves));
 for i=1:length(leaves)
     leaf=strrep(num2str(leaves{i}), ' ', '');
     [bool,p]=ismember(leaf,uIDS);
     clone_position(1,i)=p;
     label{p}=strrep(['PC',num2str(i)], ' ', '');
 end
 
 unobser=setdiff([2:1:length(uIDS)], clone_position);
 for i=1:length(unobser)
     label{unobser(i)}= strrep(['NK',num2str(i)], ' ', '');
 end
 
 uIDS1={};
 for i=1:length(uIDS)
     t=uIDS{i};
     t1=max(find(t~='0'));
     uIDS1{i} = t(1:t1);
 end

uIDS2=unique(uIDS1);


adjmatrix1=zeros(length(uIDS2),length(uIDS2));
for i=1:length(code)
    for j=1:(length(code{i})-1)
        scode1=code{i}{j};
        scode2=code{i}{j+1};
        [bool,p1]=ismember(scode1,uIDS);
        [bool,p2]=ismember(scode2,uIDS);
        scode3=uIDS1{p1};
        scode4=uIDS1{p2};
        [bool,p3]=ismember(scode3,uIDS2);
        [bool,p4]=ismember(scode4,uIDS2);
        adjmatrix1(p3,p4)=1; 
    end
end


for i=1:length(uIDS2)
    adjmatrix1(i,i)=0;
end


label1={};
clonenum={};
for i=1:length(uIDS2)
     label1{i}=strrep(['N',num2str(i)], ' ', '');
     clonenum{i}=0;
end


clone_position1=zeros(1,length(leaves));
 for i=1:length(leaves)
     leaf=strrep(num2str(leaves{i}), ' ', '');
     [bool,p]=ismember(leaf,uIDS2);
     clone_position1(1,i)=p;
     label1{p}=strrep(['PC',num2str(i)], ' ', '');
     clonenum{p}= length(clone_cells{i});
 end
 
 unobser=setdiff([1:1:length(uIDS2)], clone_position1);
 for i=1:length(unobser)
     label1{unobser(i)}= strrep(['NK',num2str(i)], ' ', '');
     clonenum{unobser(i)}=1;
 end
 
 
scode1= '0';
scode2='1';
[bool1,p1]=ismember(scode1,uIDS2);
[bool2,p2]=ismember(scode2,uIDS2);
if bool1
     adjmatrix1(p1,p2)=1;
     label1{p1}='root';
     bg=biograph(adjmatrix1,label1);
     for i=1:length(clonenum)
          set(bg.nodes(i), 'fontsize',clonenum{i});
     end
end
if ~bool1
     adjmatrix2= zeros(length(uIDS2)+1, length(uIDS2)+1);
     adjmatrix2(1,p2+1)=1;
     adjmatrix2(2:length(uIDS2)+1,2:length(uIDS2)+1)=adjmatrix1;
     label2=['root',label1];
     bg=biograph(adjmatrix2,label2);
     set(bg.nodes(1), 'fontsize',1);
     for i=1:length(clonenum)
         set(bg.nodes(i+1), 'fontsize',clonenum{i});
     end
end 
 

 set(bg.nodes,'shape','circle','color',[1,1,1],'lineColor',[0,0,1]);
 clonesize=zeros(1,length(leaves));
 for i=1:length(leaves)
     clonesize(1,i) = length(clone_cells{i});
 end

 cell_id=zeros(1,length(leaves));
 for i=1:length(leaves)
     cell_id(1,clone_cells{i})=i;
 end
 
 set(bg.nodes,'textColor',[0,0,0],'lineWidth',3);
 set(bg,'arrowSize',15,'edgeFontSize',15);
 get(bg.nodes,'position')
 view(bg);
 toc
