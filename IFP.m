%The input folder,and skip it.
DicomFolder = 'C:\Users\dell\Desktop\IFP\008_AxLAVADYNFROM3\';
filelist = dir(DicomFolder);
i=1;
    while i<=length(filelist)
        if filelist(i).isdir==1
            filelist = filelist([1:i-1 i+1:end]);   % skip folders
        else
            i=i+1;
        end
    end
    
%%% The line we got the dicom data.
line_begin = [121,147];  %The line start point. [x,y]
line_end = [127,152];    %The line end point.
 
%Divide the line into some pieces.
step_num =10; %
x_list = line_begin(1):double(abs(line_begin(1)-line_end(1)))/step_num:line_end(1);
y_list = line_begin(2):double(abs(line_begin(2)-line_end(2)))/step_num:line_end(2);

%Read the number of slice and group
info = dicominfo([DicomFolder '/' filelist(1).name]);
PositionNum = info.NumberOfTemporalPositions;
Images = info.ImagesInAcquisition;
GroupNum = Images / PositionNum; 

% The slice and group intersted
ISilce = 17;
IGroup = 6;

info = dicominfo([DicomFolder '/' filelist(ISilce + GroupNum*IGroup).name]);

%Covert the pixel into length in the real length.
ReconstructionDiameter = info.ReconstructionDiameter;
Row = info.Rows; 
Column = info.Columns;
Position = (x_list - line_begin(1))/double(Row)*double(ReconstructionDiameter);

% Read the time,and covert the time.
info = dicominfo([DicomFolder '/' filelist(ISilce).name]);
time0 = info.FileModDate;
time0 = ['01-Mar-2000 ',time0(12:19)];
para = 1.157407338420550e-05;

Time = [0];
MaxValPos = [0];
dicom = dicomread([DicomFolder '/' filelist(ISilce + GroupNum*IGroup).name]);
%%%Show the line in the dicom.
figure;imshow(dicom,[]);title('The origin image(blue line is the data we care.)');hold on;
plot(x_list,y_list,'b','LineWidth',2);
hold off;

%
display = [7;30;50;56];

figure;
for p= IGroup:1:PositionNum-1
    dicom = dicomread([DicomFolder '/' filelist(ISilce + GroupNum*p).name]);
    info = dicominfo([DicomFolder '/' filelist(ISilce + GroupNum*p).name]);
    time = info.FileModDate;

    %%%%%%%Use the bilinear interpolation to caculate the divided position.
    signal = [];
    for i =1:(step_num+1)
        x0 = floor(x_list(i)):1:(floor(x_list(i))+1); y0 = floor(y_list(i)):1:(floor(y_list(i))+1);
        z0 = [];
        for j = 1:length(x0)
            for k  = 1:length(y0)
                z0(j,k) = dicom(x0(j),y0(k));
            end
        end
        z1 = interp2(double(x0),double(y0),double(z0),x_list(i),y_list(i));
        signal  = [signal;z1];
    end
    signal = signal';
    
    %Caculte the time interval.
    time = ['01-Mar-2000 ',time(12:19)];
    Time = [Time;(datenum(time) - datenum(time0))/para];
    
    if (p==display(1)|p==display(2)|p==display(3)|p==display(4))
    %figure;
    if(p==display(1))
        subplot(2,2,1);
    end
    if(p==display(2))
        subplot(2,2,2);
    end
    if(p==display(3))
        subplot(2,2,3);
    end
    if(p==display(4))
        subplot(2,2,4);
    end
    scatter(Position,signal);xlabel('Position(mm)');ylabel('Signal intensity(a.u)');title([num2str(int8(Time(p-IGroup+1))),'s']);
    hold on; 
    end
    

    %%%%% best fits of a polynomial of the 6th degree
    p1=polyfit(Position,signal,6);
    Pos=linspace(min(Position),max(Position),100);
    fitval=polyval(p1,Pos);
    
    %plot the fitted result
    if (p==display(1)|p==display(2)|p==display(3)|p==display(4))
    plot(Pos,fitval);  
    end
    
    %Get the max value and its position.
    [val,pos] = max(fitval);
    
    %draw the max value  
    if (p==display(1)|p==display(2)|p==display(3)|p==display(4))
    stem(Pos(pos),double(val),'-.',...
       'MarkerFaceColor','red');
    end
    
    MaxValPos = [MaxValPos;Pos(pos)];
end


Time(1) = 0;
figure;scatter(Time,MaxValPos);xlabel('Time(s)');ylabel('Displacement(mm)');
%%fit the function: S(t) = S0*(1-exp(-b*t));
BETA0 = [0.5;0.01];
f_model=@(b,t)(b(1)*(1-exp(-b(2)*t)));
%opts = statset('nlinfit');
%opts.RobustWgtFun = 'bisquare';
b = nlinfit(Time,MaxValPos,f_model,BETA0);

%%Use the fitted function 
t = min(Time):0.01:max(Time);
S = b(1)*(1-exp(-b(2)*t));
S0 = b(1);
b = b(2);
hold on;plot(t,S,'r');
v0 = S0*b;
fprintf(' %s %s %s %s %s %s\n','S0:',S0,'b:',b,'v0:',v0);
