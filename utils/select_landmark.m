function [source_landmark, target_landmark] = select_landmark(source_face, source_vertex,...
                                                target_face, ...
                                                target_vertex, varargin)
% varargin: source intensity and target intensity
if size(varargin) == 0
    fig1 = figure(1); set(fig1, 'Position', [20, 100, 500, 500]);
    gpp_plot_mesh(source_face, source_vertex); title('source mesh');
    fig2 = figure(2); set(fig2, 'Position', [600, 100, 500, 500]);
    gpp_plot_mesh(target_face, target_vertex); title('target mesh')
else
    assert(length(varargin) == 2, 'only accepting single channel source and target intensity') 
    fig1 = figure(1);
    gpp_plot_mesh(source_face, source_vertex, varargin{1}); title('source mesh');
    fig2 = figure(2);
    gpp_plot_mesh(target_face, target_vertex, varargin{2}); title('target mesh')
end

stop_signal = 0;
num_landmark = 0;
source_landmark = [];
target_landmark = [];
while stop_signal == 0
    disp('In Figure 1, find the desired location then hit Return')
    rotate3d(fig1); 
    pause;
    disp('Select the landmark in Figure 1 then hit Return')
    datacursormode(fig1);
    dcm_obj1 = datacursormode(fig1);
    pause;
    info_struct1 = getCursorInfo(dcm_obj1);
    pos1 = info_struct1.Position;
    ind1 = get_index_by_location(pos1, source_vertex);
    fprintf('selected location for source %f %f %f, index %d \n', pos1(1), pos1(2), pos1(3), ind1);
    
    disp('In Figure 2, find the desired location then hit Return')
    rotate3d(fig2);
    pause;
    disp('Select the landmark in Figure 2 then hit Return')
    datacursormode(fig2);
    dcm_obj2 = datacursormode(fig2);
    pause;
    info_struct2 = getCursorInfo(dcm_obj2);
    pos2 = info_struct2.Position;
    ind2 = get_index_by_location(pos2, target_vertex); 
    fprintf('selected location for target %f %f %f, index %d \n', pos2(1), pos2(2), pos2(3), ind2);
    
    confirm_flag = input('confirm this landmark pair? yes=1 or no=0: \n');
    if confirm_flag == 1
        source_landmark = [source_landmark, ind1];
        target_landmark = [target_landmark, ind2];
        num_landmark = num_landmark + 1;
    end
    fprintf('Total number of landmark: %d', num_landmark);
    stop_signal = input('confirm finish?, yes=1 or no=0: ');
    
    
end
    

end

