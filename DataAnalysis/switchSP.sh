# 0 switch from para to ser
# 1 switch from ser to para
if [ -e  "./template.lua.ser" ]; then
mv template.lua template.lua.para
mv template.lua.ser template.lua
elif [ -e "./template.lua.para" ]; then
mv template.lua template.lua.ser
mv template.lua.para template.lua
fi
