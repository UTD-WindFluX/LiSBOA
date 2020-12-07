%08/28/2019: hor widened
function f=mkfig(varargin)

    switch nargin 
        case 0
            f=figure('Color','w','Units','normalized');

        case 1
            switch varargin{1}
                case 'max'
                    f=figure('Color','w','Units','normalized','Position', [0.0187 0.0528 0.9714 0.8602]);
                
                case 'hor'
                    f=figure('Color','w','Units','normalized','Position', [0.0055 0.3343 0.9904 0.4861]);
                case 'ver'
                    f=figure('Color','w','Units','normalized','Position',   [0.317187500000000 0.056481481481481 0.474739583333333 0.824074074074074]);
                case 'square'
                    f=figure('Color','w','Units','normalized','Position',[0.3171 0.0657 0.5534 0.8550]);
            end
    end
end