function run_color_L1CFT_ACLKS_benchmark
clc;
clear;
numColor = 10;
colorset=cell(10,1);
colorset{1}='GRAY';
colorset{2}='RGB';
colorset{3}='NRGB';
colorset{4}='HSV';
colorset{5}='OPPONENT';
colorset{6}='COPPONENT';
colorset{7}='NOPPONENT';
colorset{8}='HUE';
colorset{9}='TRGB';
colorset{10}='LAB';

numSeq = 129;
seqset=cell(numSeq, 1);
seqset{1}='Basketball';
seqset{2}='Bicycle';
seqset{3}='Biker';
seqset{4}='Bird';
seqset{5}='Board';
seqset{6}='Bolt';
seqset{7}='Boy';
seqset{8}='CarDark';
seqset{9}='CarScale';
seqset{10}='Coke';
seqset{11}='Couple';
seqset{12}='Crossing';
seqset{13}='Cup';
seqset{14}='David';
seqset{15}='David3';
seqset{16}='Deer';
seqset{17}='Diving';
seqset{18}='Doll';
seqset{19}='FaceOcc1';
seqset{20}='Football1';
seqset{21}='Girl';
seqset{22}='Girlmov';
seqset{23}='Gym';
seqset{24}='Hand';
seqset{25}='Iceskater';
seqset{26}='Ironman';
seqset{27}='Jogging1';
seqset{28}='Jogging2';
seqset{29}='Juice';
seqset{30}='Lemming';
seqset{31}='Liquor';
seqset{32}='Matrix';
seqset{33}='MotorRolling';
seqset{34}='MountainBike';
seqset{35}='Panda';
seqset{36}='Shaking';
seqset{37}='Singer1';
seqset{38}='Singer2';
seqset{39}='Skating1';
seqset{40}='Skating2';
seqset{41}='Skiing';
seqset{42}='Soccer';
seqset{43}='Subway';
seqset{44}='Sunshade';
seqset{45}='Tiger1';
seqset{46}='Tiger2';
seqset{47}='Torus';
seqset{48}='Trellis';
seqset{49}='Walking';
seqset{50}='Walking2';
seqset{51}='Woman';

seqset{52}='Airport_ce';
seqset{53}='Baby_ce';
seqset{54}='Badminton_ce1';
seqset{55}='Badminton_ce2';
seqset{56}='Basketball_ce1';
seqset{57}='Basketball_ce2';
seqset{58}='Basketball_ce3';
seqset{59}='Bike_ce1';
seqset{60}='Bike_ce2';
seqset{61}='Bikeshow_ce';
seqset{62}='Boat_ce1';
seqset{63}='Boat_ce2';
seqset{64}='Busstation_ce1';
seqset{65}='Busstation_ce2';
seqset{66}='Carchasing_ce1';
seqset{67}='Carchasing_ce3';
seqset{68}='Carchasing_ce4';
seqset{69}='Eagle_ce';
seqset{70}='Electricalbike_ce';
seqset{71}='Face_ce';
seqset{72}='Guitar_ce1';
seqset{73}='Guitar_ce2';
seqset{74}='Hurdle_ce1';
seqset{75}='Hurdle_ce2';
seqset{76}='Kite_ce1';
seqset{77}='Kite_ce2';
seqset{78}='Kite_ce3';
seqset{79}='Kobe_ce';
seqset{80}='Logo_ce';
seqset{81}='Messi_ce';
seqset{82}='Michaeljackson_ce';
seqset{83}='Motorbike_ce';
seqset{84}='Plane_ce2';
seqset{85}='Railwaystation_ce';
seqset{86}='Singer_ce1';
seqset{87}='Singer_ce2';
seqset{88}='Skating_ce1';
seqset{89}='Skating_ce2';
seqset{90}='Skiing_ce';
seqset{91}='Skyjumping_ce';
seqset{92}='Spiderman_ce';
seqset{93}='Suitcase_ce';
seqset{94}='Surf_ce1';
seqset{95}='Surf_ce2';
seqset{96}='Surf_ce3';
seqset{97}='Surf_ce4';
seqset{98}='Tennis_ce1';
seqset{99}='Tennis_ce2';
seqset{100}='Tennis_ce3';
seqset{101}='Toyplane_ce';

seqset{102}='Ball_ce1';
seqset{103}='Ball_ce2';
seqset{104}='Ball_ce3';
seqset{105}='Ball_ce4';
seqset{106}='Bee_ce';
seqset{107}='Charger_ce';
seqset{108}='Cup_ce';
seqset{109}='Face_ce2';
seqset{110}='Fish_ce1';
seqset{111}='Fish_ce2';
seqset{112}='Hand_ce1';
seqset{113}='Hand_ce2';
seqset{114}='Microphone_ce1';
seqset{115}='Microphone_ce2';
seqset{116}='Plate_ce1';
seqset{117}='Plate_ce2';
seqset{118}='Pool_ce1';
seqset{119}='Pool_ce2';
seqset{120}='Pool_ce3';
seqset{121}='Ring_ce';
seqset{122}='Sailor_ce';
seqset{123}='SuperMario_ce';
seqset{124}='TableTennis_ce';
seqset{125}='TennisBall_ce';
seqset{126}='Thunder_ce';
seqset{127}='Yo-yos_ce1';
seqset{128}='Yo-yos_ce2';
seqset{129}='Yo-yos_ce3';

str=pwd;
index_dir=strfind(str,'\');
str_temp=str(1:index_dir(end)-1);
base_path=strcat(str_temp,'\Temple-color-128');


if ~exist('result', 'dir')
    mkdir('result');
end
%default settings
%     params  = default_parameters_dat;
 for cm=2:2 %choose the color representation
    for i=1:129
        close all;
        colorModel = colorset{cm};
        fprintf('%s %d\n', colorModel, i);
        saveDir = ['result/' colorModel '/'];
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end
        %parameters according to the paper. at this point we can override
        %parameters based on the chosen kernel or feature type
        video = seqset{i};
        save_file = [saveDir video '_AKAL1CFT_' colorModel '.txt'];
        if exist(save_file,'file')
            continue;
        end
        [img_files, pos, target_sz, ground_truth, video_path, depth_path] = load_video_color(base_path, video);
        params.init_pos = pos;
        params.wsize = floor(target_sz);
        params.img_files = strcat(video_path, img_files);
        % main processing
        [positions, fps]=DCPF_main(params);
%         rects = [positions(:,2) - target_sz(2)/2, positions(:,1) - target_sz(1)/2];
%         rects(:,3) = target_sz(2);
%         rects(:,4) = target_sz(1);
        rects      = positions;
%         save_file = [saveDir video '_ScaleSSKCF_' colorModel '.txt'];
        fid = fopen(save_file, 'w+');
        num_row = size(rects, 1);
        for k = 1 : num_row
            fprintf(fid, '%.4f,%.4f,%.4f,%.4f\n', rects(k, 1), rects(k, 2), rects(k, 3), rects(k, 4));
        end
        fclose(fid);
    end
 end
 fprintf('finished...\n');