function connParams = connect_params_distrib()
    % Required general parameters
    connParams.numConnections = 1;
    connParams.toolName = 'DISCOS Update';
    % Optional general parameters
    
    % Parameters for connection 1
    connParams.conn(1).connectionName = 'DISCOS';
    connParams.conn(1).conn_host = 'discosweb.esoc.esa.int';
    % Personal access token for the DISCOSweb API
    connParams.conn(1).password = '<replace with token>';
end