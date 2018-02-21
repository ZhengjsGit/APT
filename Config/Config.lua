-- This is the manage configuration file of APT.
-- Please write down the path of the main configuration file.

-- Attention: 
--     1, Use "dofile" command to load.
--     2, Only one main configuration file can be loaded.
--     3, All the main configuration files locate at "./pkg/".
--     4, For Available cofiguration files and their instructions, 
--        see "./doc/MainConfigFile.txt".

-- Load the main configuration file
--      Please use the relative directory of manage configuration file---Config.lua
Config_Name="./template.lua"


ConfigPath=GAPS_APT_LuaConfig_ScriptPath();
CUR_DIR=string.sub(ConfigPath,0,string.find(ConfigPath,"Config.lua")-1)
dofile(CUR_DIR..Config_Name);
