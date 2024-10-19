function [idx] = GetConnectionIdx(connectionName, connParams)
% GetConnectionIdx - Gets the index for the connection matching the
%                    connection name passed in. If the connection name is
%                    not found in connParams, then an error is thrown.
%
% Usage: [idx] = GetConnectionIdx(connectionName, connParams)
%
% Inputs:
%   connectionName - The name of the connection as defined in the
%                    connParams structure.
%   connParams - (Optional) Connection parameters structure defining the
%                connections and credential files to be used. If this
%                parameter isn't entered, then it is assumed a
%                "connect_params.m" file exists.
%
% Outputs:
%   idx - Index for the connection matching the connection name.

    if nargin == 1
        connParams = connect_params;
    end
    connectionFound = false;
    for i = 1:connParams.numConnections
        if strcmp(connParams.conn(i).connectionName,connectionName)
            connectionFound = true;
            idx = i;
        end
    end
    if ~connectionFound
        error(['Could not find connection named ' connectionName]);
    end
end