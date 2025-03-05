classdef eProbe
    
    % eProbe   Realtime Processing of bioSignals
    %
    %
    % ScienceBeam
    
    properties 
        FreqResponse (1,1) logical = false % (boolean) sets monitoring of frequency response in showEEG mode 
        fs (1,1) double {mustBeMember(fs,[1000,500])}=1000 % Only takes either 500Hz or 1kHz
        duration (1,1) double = 60000 % (ms) Total runtime of program
        Channels double = 1:8 % List of channels which are going to be analysied
        timeLength (1,1) double = 5000 % (ms) Length of window
        FilterBlock = false %Window-based FIR filter, use createFilter or other MATLAB's functions to generate desired filter. If other MATLAB's functions are used, order of filter must be 180 for 1kHz and 450 for 500Hz
        
        
        
    end
    methods(Static)
        
function showEEG(obj)
    % Real time monitoring of channels.
    fs = obj.fs;
    FreqResp = obj.FreqResponse;
    duration = obj.duration;
    Ch = obj.Channels;
    timeLength = obj.timeLength;
    FilterBlock = obj.FilterBlock;
    if fs==500
        duration=duration/2;
        timeLength=timeLength/2;
    end
    if FilterBlock~=false
        timeLength = timeLength+1000 ;
    end
    
    obj.FigureProperties()
    
    ploterData = parallel.pool.DataQueue;
    ploterData.afterEach(@(x) obj.DataUpdater(FreqResp,timeLength,x,Ch,FilterBlock,fs))
    
    
    
    if fs==500
        N = floor(2*duration*fs/1000);
    else
        N = floor(duration*fs/1000);
    end
    
spmd(2) %%%% Workers
    
    
    if labindex==1
        
        NET.addAssembly('C:\RubyMind\eProbe\eBridge1.dll');  % pathway of dll
        eBridge.eBridgeDAQ.Init();
        eBridge.eBridgeDAQ.Start()
        b = zeros(100,length(Ch));
        k=1;
        w=0;
        counter=0;
        while 1
            
            readPtr=eBridge.eBridgeDAQ.GetNextReadPointer(0,1);
            b = obj.select(Ch,k,b,readPtr);
            
            if k==100
                
                labSend(b,2);
                k=0;
                counter = counter+1;
            end
            k=k+1;
            w=w+1;
            
            
            if counter == N/100
                eBridge.eBridgeDAQ.Stop();
                break;
            end
        end
    end
    
    
    
    if labindex==2
        
        antiAliasFilter5 = fir1(99,2*40/fs,'low',chebwin(100,20));
        antiAliasFilter5 = antiAliasFilter5';
        
        for r=1:length(Ch)
            dc2{r} = dsp.DCBlocker('Algorithm','FIR','Length', 100);
        end
        Data = zeros(timeLength,length(Ch)); %%%% temporary storage
        
        if fs==500
        
        DS5Data = zeros(timeLength/2,length(Ch)); %%%% DownSampled Data /5
        FinalData=zeros(timeLength/2,length(Ch));
        count=0;
        t=1;
        for cicl = 1:ceil(duration/timeLength)
            
            for a=1:timeLength/100
                pack=[];
                Receiver = labReceive(1); %Catch data from other worker
                for r=1:length(Ch)
                
                    data = dc2{r}(Receiver(:,r));
                    %%%% downsampling process %%%%
                    
                    Data(:,r) = obj.BlockFilteringDS5(data,antiAliasFilter5,a,Data(:,r),timeLength);
                    DS5Data((a-1)*50+1:a*50,r) = downsample(Data((a-1)*100+1:a*100,r),2);

                    

                    if FilterBlock ~= false
                        FinalData(:,r) = obj.BlockFiltering(FilterBlock,FinalData(:,r)...
                            ,DS5Data((a-1)*50+1:a*50,r),a,timeLength,fs);
                        cData = obj.cutData(FinalData(:,r),a,t,fs);
                    else
                        cData = DS5Data(:,r);
                        FinalData(:,r) = DS5Data(:,r);
                    end



                    if FreqResp
                        FData = cat(1,cData...
                            ,obj.powerCalculator(obj,FinalData(:,r),a,timeLength,fs));
                        pack = cat(2,pack,FData);
                        
                    else
                        pack = cat(2,pack,cData);
                    end

                end
                
                if FreqResp
                    if mod(a,2)
                        send(ploterData, pack); %Sending Data for monitoring
                    end
                else
                    send(ploterData, pack);
                end
                
                data=[]; % Clear dataRegister
                count = count+1;
                if count== duration/100
                    break
                end
                
                t=mod(t,(timeLength-1000)/100)+1;
            end
            
        end
        
        else
        
        DS5Data = zeros(timeLength/5,length(Ch)); %%%% DwonSampled Data /5
        FinalData=zeros(timeLength/5,length(Ch));
        count=0;
        t=1;
        for cicl = 1:ceil(duration/timeLength)
            
            for a=1:timeLength/100
                pack=[];
                Receiver = labReceive(1); %Catch data from other worker
                for r=1:length(Ch)
                
                    data = dc2{r}(Receiver(:,r));
                    %%%% downsampling process %%%%
                    
                    Data(:,r) = obj.BlockFilteringDS5(data,antiAliasFilter5,a,Data(:,r),timeLength);
                    DS5Data((a-1)*20+1:a*20,r) = downsample(Data((a-1)*100+1:a*100,r),5);

                    

                    if FilterBlock ~= false
                        FinalData(:,r) = obj.BlockFiltering(FilterBlock,FinalData(:,r)...
                            ,DS5Data((a-1)*20+1:a*20,r),a,timeLength,fs);
                        cData = obj.cutData(FinalData(:,r),a,t,fs);  
                    else
                        cData = DS5Data(:,r);
                        FinalData(:,r) = DS5Data(:,r);
                    end



                    if FreqResp
                        FData = cat(1,cData...
                            ,obj.powerCalculator(obj,FinalData(:,r),a,timeLength,fs));
                        pack = cat(2,pack,FData);
                        
                    else
                        pack = cat(2,pack,cData);
                    end

                end
                
                if FreqResp
                    if mod(a,2)
                        send(ploterData, pack); %Sending Data for monitoring
                    end
                else
                    send(ploterData, pack);
                end
                
                data=[]; % Clear dataRegister
                count = count+1;
                if count== duration/100
                    break
                end
                
                 t=mod(t,timeLength/100-10)+1;
            end
            
        end
        end
    end
end

end    


function corMap(obj)
    % color map of crrelations between Channels
    Ch = obj.Channels;
    timeLength = obj.timeLength;
    duration = obj.duration;
    fs = obj.fs;
    if fs==500
        duration=duration/2;
        timeLength=timeLength/2;
    end
    
    obj.FigureProperties()
    ploterData = parallel.pool.DataQueue;
    ploterData.afterEach(@(x) obj.corMapUpdater(x,timeLength,Ch,fs))
    
    
    if fs==1000
        N = floor(duration*fs/1000);
    else
        N = floor(2*duration*fs/1000);
    end
   
    
spmd(2) %%%% Workers
    
    Chnum = length(Ch);
    
    if labindex==1
        NET.addAssembly('C:\RubyMind\eProbe\eBridge1.dll');  % pathway of dll
        eBridge.eBridgeDAQ.Init();
        eBridge.eBridgeDAQ.Start()
        b = zeros(100,Chnum);
        k=1;
        w=0;
        counter=0;
        while 1
            
            readPtr=eBridge.eBridgeDAQ.GetNextReadPointer(0,1);
            b = obj.select(Ch,k,b,readPtr);
            
            if k==100
                
                labSend(b,2);
                k=0;
                counter = counter+1;
            end
            k=k+1;
            w=w+1;
            
            
            if counter == N/100
                eBridge.eBridgeDAQ.Stop();
                break;
            end
        end
    end
    
    
    
    if labindex==2
        
        antiAliasFilter5 = fir1(99,2*45/fs,'low',chebwin(100,50));
        antiAliasFilter5 = antiAliasFilter5';
        
        for r=1:Chnum
            dc2{r} = dsp.DCBlocker('Algorithm','FIR','Length',100);
        end
        Data = zeros(timeLength,Chnum); %%%% temporary storage
        
        if fs==500
        
        argData = zeros(timeLength/2,Chnum);
        DS5Data = zeros(timeLength/2,Chnum); %%%% DwonSampled Data /5
        Dcov = zeros(Chnum);
        count=0;
        
        for cicl = 1:ceil(duration/timeLength)
            
            for a=1:timeLength/100
                pack=[];
                Receiver = labReceive(1); %Catch data from other worker
                for r=1:Chnum
                
                    data = dc2{r}(Receiver(:,r));
                    %%%% downsampling process %%%%
                    
                    Data(:,r) = obj.BlockFilteringDS5(data,antiAliasFilter5,a,Data(:,r),timeLength);
                    DS5Data((a-1)*50+1:a*50,r) = downsample(Data((a-1)*100+1:a*100,r),2);
                    argData(:,r) = obj.arngData(DS5Data(:,r),a,timeLength,fs);
                    
                end
                
                for i=1:Chnum
                    for j=1:Chnum
                        Dcov(i,j) = obj.covCalculator(argData(:,i),argData(:,j));
                    end
                end
                
                
                send(ploterData, Dcov);
                
                data=[]; % Clear dataRegister
                count = count+1;
                if count== duration/100
                    break
                end
                
                 
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
        
        argData = zeros(timeLength/5,Chnum);
        DS5Data = zeros(timeLength/5,Chnum); %%%% DwonSampled Data /5
        Dcov = zeros(Chnum);
        count=0;
        
        for cicl = 1:ceil(duration/timeLength)
            
            for a=1:timeLength/100
                pack=[];
                Receiver = labReceive(1); %Catch data from other worker
                for r=1:Chnum
                
                    data = dc2{r}(Receiver(:,r));
                    %%%% downsampling process %%%%
                    
                    Data(:,r) = obj.BlockFilteringDS5(data,antiAliasFilter5,a,Data(:,r),timeLength);
                    DS5Data((a-1)*20+1:a*20,r) = downsample(Data((a-1)*100+1:a*100,r),5);
                    argData(:,r) = obj.arngData(DS5Data(:,r),a,timeLength,fs);
                    
                end
                
                for i=1:Chnum
                    for j=1:Chnum
                        Dcov(i,j) = obj.covCalculator(argData(:,i),argData(:,j));
                    end
                end
                
                
                send(ploterData, Dcov);
                
                data=[]; % Clear dataRegister
                count = count+1;
                if count== duration/100
                    break
                end
                
                 
            end
            
        end
        end
    end
end


end    


function chCohere(ch1,ch2,obj)
    % Monitoring of Coherence, cross-correlation, and frequency response of two channels
    timeLength = obj.timeLength;
    duration = obj.duration;
    fs = obj.fs;
    if fs==500
        duration=duration/2;
        timeLength=timeLength/2;
    end
    
    ploterData = parallel.pool.DataQueue;
    ploterData.afterEach(@(x) obj.cohDataUpdater(x,timeLength,fs))
    
    if fs==1000
        N = floor(duration*fs/1000);
    else
        N = floor(2*duration*fs/1000);
    end
    
spmd(2) %%%% Workers
    
    
    
    if labindex==1
        
        NET.addAssembly('C:\RubyMind\eProbe\eBridge1.dll');  % pathway of dll
        eBridge.eBridgeDAQ.Init();
        eBridge.eBridgeDAQ.Start()
        b = zeros(100,2);
        k=1;
        w=0;
        counter=0;
        while 1
            
            
            readPtr=eBridge.eBridgeDAQ.GetNextReadPointer(0,1);
            b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch1-1,readPtr);
            b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch2-1,readPtr);
            
            
            if k==100
                
                labSend(b,2);
                k=0;
                counter = counter+1;
            end
            k=k+1;
            w=w+1;
            
            
            if counter == N/100
                eBridge.eBridgeDAQ.Stop();
                break;
            end
        end
    end
    
    
    
    if labindex==2
        
        antiAliasFilter5 = fir1(99,2*45/fs,'low',chebwin(100,50));
        antiAliasFilter5 = antiAliasFilter5';
        
        for r=1:2
            dc2{r} = dsp.DCBlocker('Algorithm','FIR','Length', 100);
        end
        Data = zeros(timeLength,2); %%%% temporary storage
        
        if fs==500
        
        DS5Data = zeros(timeLength/2,2); %%%% DwonSampled Data /5
        
        count=0;
        
        for cicl = 1:ceil(duration/timeLength)
            
            for a=1:timeLength/100
                pack=[];
                Receiver = labReceive(1); %Catch data from other worker
                for r=1:2
                
                    data = dc2{r}(Receiver(:,r));
                    %%%% downsampling process %%%%
                    
                    Data(:,r) = obj.BlockFilteringDS5(data,antiAliasFilter5,a,Data(:,r),timeLength);
                    DS5Data((a-1)*50+1:a*50,r) = downsample(Data((a-1)*100+1:a*100,r),2);
                    
                end
                
                arngData1 = obj.arngData(DS5Data(:,1),a,timeLength,fs);
                arngData2 = obj.arngData(DS5Data(:,2),a,timeLength,fs);
                FData1 = obj.fftCalculator(arngData1,timeLength,fs);
                FData2 = obj.fftCalculator(arngData2,timeLength,fs);
                Fcoh12 = obj.coherence_MVDR(arngData1,arngData2,50,1000);
                Fcor12 = xcorr(arngData1,arngData2);
                Dcov = obj.covCalculator(arngData1,arngData2);
                pack = cat(1,FData1,FData2,Fcoh12,Fcor12,Dcov);
                send(ploterData, pack);
                
                data=[]; % Clear dataRegister
                count = count+1;
                if count== duration/100
                    break
                end
                
                 
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
        
        DS5Data = zeros(timeLength/5,2); %%%% DwonSampled Data /5
        
        count=0;
        
        for cicl = 1:ceil(duration/timeLength)
            
            for a=1:timeLength/100
                pack=[];
                Receiver = labReceive(1); %Catch data from other worker
                for r=1:2
                
                    data = dc2{r}(Receiver(:,r));
                    %%%% downsampling process %%%%
                    
                    Data(:,r) = obj.BlockFilteringDS5(data,antiAliasFilter5,a,Data(:,r),timeLength);
                    DS5Data((a-1)*20+1:a*20,r) = downsample(Data((a-1)*100+1:a*100,r),5);
                    
                end
                
                arngData1 = obj.arngData(DS5Data(:,1),a,timeLength,fs);
                arngData2 = obj.arngData(DS5Data(:,2),a,timeLength,fs);
                FData1 = obj.fftCalculator(arngData1,timeLength,fs);
                FData2 = obj.fftCalculator(arngData2,timeLength,fs);
                Fcoh12 = obj.coherence_MVDR(arngData1,arngData2,50,1000);
                Fcor12 = xcorr(arngData1,arngData2);
                
                Dcov = obj.covCalculator(arngData1,arngData2);
                pack = cat(1,FData1,FData2,Fcoh12,Fcor12,Dcov);
                send(ploterData, pack);
                
                data=[]; % Clear dataRegister
                count = count+1;
                if count== duration/100
                    break
                end
                
                 
            end
            
        end
        end
    end
end


end  


function myFilter = createFilter(myband,type,obj)
    % Inputs: myband {'alpha', 'beta', 'delta', 'theta', a, [a,b]} type {[], 'bandpass', 'stop', 'lowpass', 'highpass'}
    fs = obj.fs;
    if length(myband)==1
        if isempty(type)
            error('please enter type of your filter');
        else
            band=myband;
        end
        
    elseif strcmp(type,'bandpass')
        band = myband;
        
    elseif strcmp(type,'stop')
        band = myband;
    
    elseif strcmp(myband,'beta')
        band = [13,32]; type='bandpass';
    
    elseif strcmp(myband,'alpha')
        band=[8,13]; type='bandpass';
    
    elseif strcmp(myband,'theta')
        band=[4,8]; type='bandpass';
    
    elseif strcmp(myband,'delta')
        band=4; type='lowpass';
                
    else
        
        error('Invalid syntax');
        
    end
    
    if fs==1000
        Fs = 200;
        N =180;
    elseif fs==500
        Fs = 250;
        N =450;
    else
        error('Choose either 500Hz or 1kHz')
    end
    
    for i=1:length(band)
        if band(i)>= Fs/2
            error('Band is greater than Nyquist frequency');
        end
    end
    if strcmp(type,'lowpass')
        myFilter = fir1(N-1,band/Fs*2,'low',chebwin(N,40));
    elseif strcmp(type,'highpass')
        
        myFilter = firpm(N,[0 band 1.1*band Fs/2]/(Fs/2),[0 0 1 1],[40 1]);
        myFilter=myFilter(1:end-1);

    elseif strcmp(type,'bandpass')
        myFilter = fir1(N-1,band/Fs*2,'bandpass',chebwin(N,40));
    elseif strcmp(type,'stop')
        if fs==1000
            band = [0.92*band(1),band,1/.92*band(2)];
            myFilter = firpm(N-2,[0 band  Fs/2]/(Fs/2),[1 1 0 0 1 1],[1 100 1]);
        elseif fs==500
            band = [0.96*band(1),band,1/.96*band(2)];
            myFilter = firpm(N-2,[0 band  Fs/2]/(Fs/2),[1 1 0 0 1 1],[1 100 1]);
        end
        myFilter=[myFilter,0];
    else
        error('Invalid input');
           
    end
    myFilter = myFilter';
    
    
end


    end
    
    methods (Hidden=true, Static=true, Sealed=true)
function FigureProperties()
    
    figure('Name','EEG-matProb','NumberTitle','off','position',get(0,'screensize'));
    
end


function BandPower = powerCalculator(obj,Data,a,timeLength,fs)
    
    x = obj.arngData(Data,a,timeLength,fs);
    
    X = fft(x,1024);
    if fs==1000
        X = X(1:256); %%%% Until 50Hz
    else
        X = X(1:205);
    end
    
    if fs==500
        Px = abs(X).^2*2/fs/length(x)/2;
    else
        Px = abs(X).^2*2/fs/length(x)/5;
    end

    N1=floor(0.4/50*length(X)); % Band Selection
    N2=floor(4/50*length(X));
    N3=floor(8/50*length(X));
    N4=floor(13/50*length(X));
    N5=floor(32/50*length(X));
    d=sum(Px(N1:N2)); %delta
    t=sum(Px(N2:N3)); %theta
    a=sum(Px(N3:N4)); %alpha
    b=sum(Px(N4:N5)); %beta
    
    BandPower = [d,t,a,b]';

end


function DataUpdater(FreqResp,timeLength,q,Ch,FilterBank,fs)
    
    if FilterBank==false
        
       L = timeLength;
       
    else
        L = timeLength-1000;
    end
    if fs==1000
        h=5;
    elseif fs==500
        L=L*2;
        h=4;
    end
    
    L=L/h;
    Chnum = length(Ch);
    bias = zeros(L,Chnum);
    drawnow nocallbacks
    for i=1:length(Ch)
        bias(:,i) = -2000*(Chnum-1)+(i-1)*4000;
    end
    for i=1:Chnum
        k{i} = strcat('Ch ',num2str(Ch(i)));
    end
    if FreqResp
        
        subplot(2,5,[1,2,3,4,6,7,8,9]);
        
        
        plot(q(1:L,:)+bias)
        set(gca,'XTick',[],'YTick',4000*((1:Chnum)-1)-2000*(Chnum-1),...
                'Ylim',[-2000*Chnum,2000*Chnum],'Yticklabels',k,'Xlim',[1,L],...
                'XTick',L/2,'Xticklabels',['time window ',num2str(h*L),'ms']);
        
        subplot(2,5,[5,10])
        barh(1:Chnum,q(1+L:end,:))
        set(gca,'YTick',[],'Ylim',[0.5,Chnum+.5]);
        title(['\fontsize{16} {\color[rgb]{0 0.4470 0.7410}{\delta} '...
            '\color[rgb]{0.8500 0.3250 0.0980}{\theta} \color[rgb]{0.9290 0.6940 0.1250}'...
            '{\alpha} \color[rgb]{0.4940 0.1840 0.5560}{\beta}} '])
        
    else
            
        plot(q+bias);

        set(gca,'XTick',[],'YTick',4000*((1:Chnum)-1)-2000*(Chnum-1),...
            'Ylim',[-2000*Chnum,2000*Chnum],'Yticklabels',k,'Xlim',[1,L],...
            'XTick',L/2,'Xticklabels',['time window ',num2str(L*h),'ms']);
        
    end
    
        
        
end


function FilterData = BlockFilteringDS5(data,antiAliasFilter5,a,Data,timeLength)

    batchFilter = conv(antiAliasFilter5,data);
    if a==timeLength/100

        Data(100*a-1:a*100)=0;
        Data((a-1)*100+1:a*100) =Data((a-1)*100+1:a*100)+ batchFilter(1:100);
        Data(1:99) = batchFilter(101:199);
    else
        Data(a*100:a*100+99)=0;
        Data((a-1)*100+1:a*100+99) =Data((a-1)*100+1:a*100+99)+ batchFilter;
    end
    FilterData = Data;

end


function FilterData = BlockFiltering(myFilter,fData,data,a,timeLength,fs)
    
    if fs==1000
        batchFilter = conv(myFilter,data);
        bsize = length(batchFilter);
        if (a-1)*20+bsize>timeLength/5
            l = mod((a-1)*20+bsize,timeLength/5);
            if l<20
                dl = 20-l;
                fData(end-dl+1:end)=0;fData(1:l)=0;
            else
                fData(l-20+1:l) = 0;
            end
            fData(20*(a-1)+1:end) = fData(20*(a-1)+1:end) + batchFilter(1:end-l);
            fData(1:l) = fData(1:l) + batchFilter(end-l+1:end);

        else
            fData((a-2)*20+bsize+1:(a-1)*20+bsize) = 0;
            fData((a-1)*20+1:(a-1)*20+bsize) = fData((a-1)*20+1:(a-1)*20+bsize)...
                +batchFilter;
        end
    
    elseif fs==500
        
        batchFilter = conv(myFilter,data);
        bsize = length(batchFilter);
        if (a-1)*50+bsize>timeLength/2
            l = mod((a-1)*50+bsize,timeLength/2);
            if l<50
                dl = 50-l;
                fData(end-dl+1:end)=0;fData(1:l)=0;
            else
                fData(l-50+1:l) = 0;
            end
            fData(50*(a-1)+1:end) = fData(50*(a-1)+1:end) + batchFilter(1:end-l);
            fData(1:l) = fData(1:l) + batchFilter(end-l+1:end);

        else
            fData((a-2)*50+bsize+1:(a-1)*50+bsize) = 0;
            fData((a-1)*50+1:(a-1)*50+bsize) = fData((a-1)*50+1:(a-1)*50+bsize)...
                +batchFilter;
        end
        
    end
    FilterData = fData;
end
        

function adpData = cutData(Data,a,t,fs)
    if fs==1000
        L=length(Data)/20;
        x = cat(1,Data(20*a:end),Data(1:a*20-1));
        x = x(200+1:end);
        adpData = cat(1,x((L-10-t+1)*20+1:end),x(1:(L-10-t+1)*20));
    else
        L=length(Data)/50;
        x = cat(1,Data(50*a:end),Data(1:a*50-1));
        x = x(500+1:end);
        adpData = cat(1,x((L-10-t+1)*50+1:end),x(1:(L-10-t+1)*50));
    end
end


function Dcov = covCalculator(data1,data2)
    sigcor = cov(data1,data2);
    Dcov = abs(sigcor(1,2))/sqrt(sigcor(1,1)*sigcor(2,2));
end


function cohDataUpdater(x,timeLength,fs)

    if fs==500
        k=2;
        l = floor(timeLength/5);
    else
        k=5;
        l = floor(timeLength/4);
    end
    
    FData1 = x(1:l);
    FData2 = x(l+1:2*l);
    Fcor12 = x(l*2+250+1:end-1);
    f = linspace(0,50,length(FData1));
    if fs==1000
        Fcoh12 = x(l*2+1:l*2+250);
        ff = linspace(0,50,250);
    else
        Fcoh12 = x(l*2+1:l*2+200);
        ff = linspace(0,50,200);
    end
    lag = (-timeLength/k+1:timeLength/k-1)*k;
    drawnow nocallbacks
    drawnow nocallbacks 
    subplot(2,2,[1,2])
    plot(f,[FData1,FData2])
    ylabel('Freq Response')
    xlabel('Hz')
    xlim([0,50])
    
    subplot(2,2,3)
    plot(ff,Fcoh12)
    axis([0,50,0,1])
    ylabel('Coherence')
    xlabel('Hz')
    subplot(2,2,4)
    plot(lag,Fcor12)
    xlabel('lag(ms)');
    ylabel('Correlation');
    title(['cor= ',num2str(x(end))])
end


function adjData = arngData(Data,a,timeLength,fs)

    if fs==500
        k=2;
    else
        k=5;
    end
    if a==timeLength/100
        adjData = cat(1,Data(100/k:end),Data(1:100/k-1));
    else
        adjData = cat(1,Data((a+1)*100/k:end),Data(1:(a+1)*100/k-1));
    end
end


function MSC =coherence_MVDR(x1,x2,L,K)

% This program computes the coherence function between 2 signals 
% x1 and x2 with the MVDR method.

% L is the length of MVDR filter or window length
% K is the resolution (the higher K, the better the resolution)

%initialization
    xx1     = zeros(L,1);
    xx2     = zeros(L,1);
    r11     = zeros(L,1);
    r22     = zeros(L,1);
    r12     = zeros(L,1);
    r21     = zeros(L,1);

    %construction of the Fourier Matrix
    F       = zeros(L,K);
    l       = [0:L-1]';
    f       = exp(2*pi*l*1i/K);
    for k = 0:K-1
        F(:,k+1) = f.^k;
    end
    F       = F/sqrt(L);

    %number of samples, equal to the lenght of x1 and x2
    n       = length(x1);

    for i = 1:n
        xx1 = [x1(i);xx1(1:L-1)];
        xx2 = [x2(i);xx2(1:L-1)];
        r11 = r11 + xx1*conj(xx1(1));
        r22 = r22 + xx2*conj(xx2(1));
        r12 = r12 + xx1*conj(xx2(1));
        r21 = r21 + xx2*conj(xx1(1));
    end
    r11 = r11/n;
    r22 = r22/n;
    r12 = r12/n;
    r21 = r21/n;
    %
    R11 = toeplitz(r11);
    R22 = toeplitz(r22);
    R12 = toeplitz(r12,conj(r21));
    %
    %for regularization
    Dt1     = 0.01*r11(1)*diag(diag(ones(L)));
    Dt2     = 0.01*r22(1)*diag(diag(ones(L)));
    %
    Ri11    = inv(R11 + Dt1);
    Ri22    = inv(R22 + Dt2);
    Rn12    = Ri11*R12*Ri22;
    %
    Si11    = real(diag(F'*Ri11*F));
    Si22    = real(diag(F'*Ri22*F));
    S12     = diag(F'*Rn12*F);
    %
    %Magnitude squared coherence function
    MSC     = real(S12.*conj(S12))./(Si11.*Si22);
    MSC = MSC(1:end/4);
end


function FData = fftCalculator(data,timeLength,fs)
    
    if fs==500
        k=2;
    else
        k=5;
    end
    Fdata = abs(fft(data,timeLength))/timeLength*k;
    if fs==500
        FData = Fdata(1:floor(timeLength/5));
    else
        FData = Fdata(1:floor(timeLength/4));
    end
end


function corMapUpdater(x,timeLength,Ch,fs)
    if fs==500
        timeLength=timeLength*2;
    end
    Chnum = length(Ch);
    for i=1:Chnum
        k{i} = strcat('Ch ',num2str(Ch(i)));
    end
    drawnow nocallbacks
    
    image(x*256);
    title(['Correlation colormap time window ',num2str(timeLength),'ms'])
    set(gca,'YTick',1:Chnum,'Yticklabels',k,'XTick',1:Chnum,'Xticklabels',k);
end


function B = select(ch,k,b,readPtr)
    

    ch=ch-1;
    if length(ch)==1
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
    elseif length(ch)==2
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
    elseif length(ch)==3
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
    elseif length(ch)==4
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
    elseif length(ch)==5
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
    elseif length(ch)==6
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
    elseif length(ch)==7
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
    elseif length(ch)==8
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
    elseif length(ch)==9
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
    elseif length(ch)==10
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
    elseif length(ch)==11
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
    elseif length(ch)==12
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
    elseif length(ch)==13
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
    elseif length(ch)==14
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
    elseif length(ch)==15
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
    elseif length(ch)==16
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
    elseif length(ch)==17
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
    elseif length(ch)==18
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
    elseif length(ch)==19
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
    elseif length(ch)==20
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
    elseif length(ch)==21
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
    elseif length(ch)==22
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
        b(k,22)=eBridge.eBridgeDAQ.GetScaledData(ch(22),readPtr);
    elseif length(ch)==23
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
        b(k,22)=eBridge.eBridgeDAQ.GetScaledData(ch(22),readPtr);
        b(k,23)=eBridge.eBridgeDAQ.GetScaledData(ch(23),readPtr);
    elseif length(ch)==24
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
        b(k,22)=eBridge.eBridgeDAQ.GetScaledData(ch(22),readPtr);
        b(k,23)=eBridge.eBridgeDAQ.GetScaledData(ch(23),readPtr);
        b(k,24)=eBridge.eBridgeDAQ.GetScaledData(ch(24),readPtr);
    elseif length(ch)==25
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
        b(k,22)=eBridge.eBridgeDAQ.GetScaledData(ch(22),readPtr);
        b(k,23)=eBridge.eBridgeDAQ.GetScaledData(ch(23),readPtr);
        b(k,24)=eBridge.eBridgeDAQ.GetScaledData(ch(24),readPtr);
        b(k,25)=eBridge.eBridgeDAQ.GetScaledData(ch(25),readPtr);
    elseif length(ch)==26
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
        b(k,22)=eBridge.eBridgeDAQ.GetScaledData(ch(22),readPtr);
        b(k,23)=eBridge.eBridgeDAQ.GetScaledData(ch(23),readPtr);
        b(k,24)=eBridge.eBridgeDAQ.GetScaledData(ch(24),readPtr);
        b(k,25)=eBridge.eBridgeDAQ.GetScaledData(ch(25),readPtr);
        b(k,26)=eBridge.eBridgeDAQ.GetScaledData(ch(26),readPtr);
    elseif length(ch)==27
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
        b(k,22)=eBridge.eBridgeDAQ.GetScaledData(ch(22),readPtr);
        b(k,23)=eBridge.eBridgeDAQ.GetScaledData(ch(23),readPtr);
        b(k,24)=eBridge.eBridgeDAQ.GetScaledData(ch(24),readPtr);
        b(k,25)=eBridge.eBridgeDAQ.GetScaledData(ch(25),readPtr);
        b(k,26)=eBridge.eBridgeDAQ.GetScaledData(ch(26),readPtr);
        b(k,27)=eBridge.eBridgeDAQ.GetScaledData(ch(27),readPtr);
    elseif length(ch)==28
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
        b(k,22)=eBridge.eBridgeDAQ.GetScaledData(ch(22),readPtr);
        b(k,23)=eBridge.eBridgeDAQ.GetScaledData(ch(23),readPtr);
        b(k,24)=eBridge.eBridgeDAQ.GetScaledData(ch(24),readPtr);
        b(k,25)=eBridge.eBridgeDAQ.GetScaledData(ch(25),readPtr);
        b(k,26)=eBridge.eBridgeDAQ.GetScaledData(ch(26),readPtr);
        b(k,27)=eBridge.eBridgeDAQ.GetScaledData(ch(27),readPtr);
        b(k,28)=eBridge.eBridgeDAQ.GetScaledData(ch(28),readPtr);
    elseif length(ch)==29
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
        b(k,22)=eBridge.eBridgeDAQ.GetScaledData(ch(22),readPtr);
        b(k,23)=eBridge.eBridgeDAQ.GetScaledData(ch(23),readPtr);
        b(k,24)=eBridge.eBridgeDAQ.GetScaledData(ch(24),readPtr);
        b(k,25)=eBridge.eBridgeDAQ.GetScaledData(ch(25),readPtr);
        b(k,26)=eBridge.eBridgeDAQ.GetScaledData(ch(26),readPtr);
        b(k,27)=eBridge.eBridgeDAQ.GetScaledData(ch(27),readPtr);
        b(k,28)=eBridge.eBridgeDAQ.GetScaledData(ch(28),readPtr);
        b(k,29)=eBridge.eBridgeDAQ.GetScaledData(ch(29),readPtr);
    elseif length(ch)==30
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
        b(k,22)=eBridge.eBridgeDAQ.GetScaledData(ch(22),readPtr);
        b(k,23)=eBridge.eBridgeDAQ.GetScaledData(ch(23),readPtr);
        b(k,24)=eBridge.eBridgeDAQ.GetScaledData(ch(24),readPtr);
        b(k,25)=eBridge.eBridgeDAQ.GetScaledData(ch(25),readPtr);
        b(k,26)=eBridge.eBridgeDAQ.GetScaledData(ch(26),readPtr);
        b(k,27)=eBridge.eBridgeDAQ.GetScaledData(ch(27),readPtr);
        b(k,28)=eBridge.eBridgeDAQ.GetScaledData(ch(28),readPtr);
        b(k,29)=eBridge.eBridgeDAQ.GetScaledData(ch(29),readPtr);
        b(k,30)=eBridge.eBridgeDAQ.GetScaledData(ch(30),readPtr);
    elseif length(ch)==31
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
        b(k,22)=eBridge.eBridgeDAQ.GetScaledData(ch(22),readPtr);
        b(k,23)=eBridge.eBridgeDAQ.GetScaledData(ch(23),readPtr);
        b(k,24)=eBridge.eBridgeDAQ.GetScaledData(ch(24),readPtr);
        b(k,25)=eBridge.eBridgeDAQ.GetScaledData(ch(25),readPtr);
        b(k,26)=eBridge.eBridgeDAQ.GetScaledData(ch(26),readPtr);
        b(k,27)=eBridge.eBridgeDAQ.GetScaledData(ch(27),readPtr);
        b(k,28)=eBridge.eBridgeDAQ.GetScaledData(ch(28),readPtr);
        b(k,29)=eBridge.eBridgeDAQ.GetScaledData(ch(29),readPtr);
        b(k,30)=eBridge.eBridgeDAQ.GetScaledData(ch(30),readPtr);
        b(k,31)=eBridge.eBridgeDAQ.GetScaledData(ch(31),readPtr);
    elseif length(ch)==32
        b(k,1)=eBridge.eBridgeDAQ.GetScaledData(ch(1),readPtr);
        b(k,2)=eBridge.eBridgeDAQ.GetScaledData(ch(2),readPtr);
        b(k,3)=eBridge.eBridgeDAQ.GetScaledData(ch(3),readPtr);
        b(k,4)=eBridge.eBridgeDAQ.GetScaledData(ch(4),readPtr);
        b(k,5)=eBridge.eBridgeDAQ.GetScaledData(ch(5),readPtr);
        b(k,6)=eBridge.eBridgeDAQ.GetScaledData(ch(6),readPtr);
        b(k,7)=eBridge.eBridgeDAQ.GetScaledData(ch(7),readPtr);
        b(k,8)=eBridge.eBridgeDAQ.GetScaledData(ch(8),readPtr);
        b(k,9)=eBridge.eBridgeDAQ.GetScaledData(ch(9),readPtr);
        b(k,10)=eBridge.eBridgeDAQ.GetScaledData(ch(10),readPtr);
        b(k,11)=eBridge.eBridgeDAQ.GetScaledData(ch(11),readPtr);
        b(k,12)=eBridge.eBridgeDAQ.GetScaledData(ch(12),readPtr);
        b(k,13)=eBridge.eBridgeDAQ.GetScaledData(ch(13),readPtr);
        b(k,14)=eBridge.eBridgeDAQ.GetScaledData(ch(14),readPtr);
        b(k,15)=eBridge.eBridgeDAQ.GetScaledData(ch(15),readPtr);
        b(k,16)=eBridge.eBridgeDAQ.GetScaledData(ch(16),readPtr);
        b(k,17)=eBridge.eBridgeDAQ.GetScaledData(ch(17),readPtr);
        b(k,18)=eBridge.eBridgeDAQ.GetScaledData(ch(18),readPtr);
        b(k,19)=eBridge.eBridgeDAQ.GetScaledData(ch(19),readPtr);
        b(k,20)=eBridge.eBridgeDAQ.GetScaledData(ch(20),readPtr);
        b(k,21)=eBridge.eBridgeDAQ.GetScaledData(ch(21),readPtr);
        b(k,22)=eBridge.eBridgeDAQ.GetScaledData(ch(22),readPtr);
        b(k,23)=eBridge.eBridgeDAQ.GetScaledData(ch(23),readPtr);
        b(k,24)=eBridge.eBridgeDAQ.GetScaledData(ch(24),readPtr);
        b(k,25)=eBridge.eBridgeDAQ.GetScaledData(ch(25),readPtr);
        b(k,26)=eBridge.eBridgeDAQ.GetScaledData(ch(26),readPtr);
        b(k,27)=eBridge.eBridgeDAQ.GetScaledData(ch(27),readPtr);
        b(k,28)=eBridge.eBridgeDAQ.GetScaledData(ch(28),readPtr);
        b(k,29)=eBridge.eBridgeDAQ.GetScaledData(ch(29),readPtr);
        b(k,30)=eBridge.eBridgeDAQ.GetScaledData(ch(30),readPtr);
        b(k,31)=eBridge.eBridgeDAQ.GetScaledData(ch(31),readPtr);
        b(k,32)=eBridge.eBridgeDAQ.GetScaledData(ch(32),readPtr);
    else
        
        error('a')
    end 
    B=b;
end

    end
end

