
st = 1;
t=tcpip('horizons.jpl.nasa.gov',6775);
buffSize = 1000000;
set(t,'InputBufferSize',buffSize);
set(t,'OutputBufferSize',buffSize);
fopen(t);
pause(3);
fid = fopen('horizons_sample.txt');
while(~feof(fid))
    if(t.BytesAvailable > 0)
        text = fread(t,[1,t.BytesAvailable]);
        native2unicode(text)
        line = fgetl(fid)
        fprintf(t, line);
        disp('==============================');
    end
end
% fclose(t);

% 
% fprintf(t, 'major-bodies');
% % pause(st);
% fprintf(t, '');
% % pause(st);
% fprintf(t, 'SATURN');
% % pause(st);
% fprintf(t, '699');
% % pause(st);
% fprintf(t, 'E');
% % pause(st);
% fprintf(t, 'v');
% % pause(st);
% fprintf(t, '500@10');
% % pause(st);
% fprintf(t, 'y');
% % pause(st);
% fprintf(t, 'eclip');
% % pause(st);
% fprintf(t, '2019-Nov-07 00:00');
% % pause(st);
% fprintf(t, '2019-Nov-07 23:59');
% % pause(st);
% fprintf(t, '10m');
% % pause(st);
% fprintf(t, 'n');
% % pause(st);
% fprintf(t, 'J2000');
% % pause(st);
% fprintf(t, '1');
% % pause(st);
% fprintf(t, '2');
% % pause(st);
% fprintf(t, 'YES');
% % pause(st);
% fprintf(t, 'YES');
% % pause(st);
% fprintf(t, '3');
% % pause(st);
% fprintf(t, 'M');
% % pause(st);
% fprintf(t, 'mjz52@cornell.edu');
% % pause(st);
% fprintf(t, 'yes');
% % pause(st);
% fprintf(t, 'q');
% % pause(st);
% fclose(t);
% 
% 






%file = 'test_HORIZONS.run';
% 
% 
% 
% [status, results] = dos('telnet horizons.jpl.nasa.gov 6775 &');
% dos('major-bodies');
