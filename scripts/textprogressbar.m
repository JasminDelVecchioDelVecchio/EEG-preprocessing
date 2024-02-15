classdef textprogressbar < handle
    properties
        strCR
        strOut
        strPercentageLength
        strDotsMaximum
        min
        max
    end
    methods
        
        function obj = textprogressbar(text)
            % Vizualization parameters
            obj.strPercentageLength = 10;   %   Length of percentage string (must be >5)
            obj.strDotsMaximum      = 10;   %   The total number of dots in a progress bar
            obj.min = 0;
            obj.max = 100;
            fprintf('%s',text);
            obj.strCR = -1;
        end
        
        function updatebar(obj,pctg)
            if isnumeric(pctg)
                pctg = pctg/100*(obj.max-obj.min) + obj.min;
                pctg = floor(pctg);
                percentageOut = [num2str(pctg) '%%'];
                percentageOut = [percentageOut repmat(' ',1,obj.strPercentageLength-length(percentageOut)-1)];
                nDots = floor(pctg/100*obj.strDotsMaximum);
                dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,obj.strDotsMaximum-nDots) ']'];
                obj.strOut = [percentageOut dotOut];
                % Print it on the screen
                if obj.strCR == -1,
                    % Don't do carriage return during first run
                    fprintf(obj.strOut);
                else
                    % Do it during all the other runs
                    fprintf([obj.strCR obj.strOut]);
                end
                
                % Update carriage return
                obj.strCR = repmat('\b',1,length(obj.strOut)-1);
                
            else
                obj.strOut = [pctg];
                
                % Print it on the screen
                if obj.strCR == -1,
                    % Don't do carriage return during first run
                    fprintf(obj.strOut);
                else
                    % Do it during all the other runs
                    fprintf([obj.strCR obj.strOut]);
                end
                
                % Update carriage return
                obj.strCR = repmat('\b',1,length(obj.strOut));

            end
                        
        end
        
        function endprogressbar(obj,text)
            obj.strCR = -1;
            fprintf([text '\n']);
        end
    end
end