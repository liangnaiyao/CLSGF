function data=data_norm(data,data_g)
if data_g==1
    num=size(data{1},2);
    v=length(data);
    for i = 1 :v  
        for  j = 1:num 
            data{i}(:,j) = ( data{i}(:,j)- mean( data{i}(:,j) ) ) / std( data{i}(:,j) ) ; 
        end 
    end
elseif data_g==2   
   for v = 1:length(data)
     A = mapminmax(data{v},0,1); 
     data{v} = A; 
   end    
elseif data_g==3
    for v = 1:length(data)
     data{v} = data{v}*diag(sparse(1./sqrt(sum(data{v}.^2))));
    end  
 end    