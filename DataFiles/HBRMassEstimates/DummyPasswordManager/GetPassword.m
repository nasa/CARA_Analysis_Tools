function [password] = GetPassword(connectionName, ~, connParams)
    connidx = GetConnectionIdx(connectionName, connParams);
    password = connParams.conn(connidx).password;
end

