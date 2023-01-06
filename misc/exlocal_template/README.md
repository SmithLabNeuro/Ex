# Ex_local
Smith Lab Rig Computer Ex_local Default Template<BR>

This is an example "Ex_local" directory. You can move this directory into the same level as the Ex directory and rename it as "Ex_local". So, it should be configured like this:<BR>

/home/username/Ex<BR>
/home/username/Ex_local<BR>

This command should do what is necessary:<BR>

cd /home/username/Ex<BR>
cp -r misc/exlocal_template ../Ex_local<BR>

Ex will automatically add this to your path and you can place your experiment files and XML configuration files here. This enables you to keep the "Ex" directory synchronized from github while not overwriting your own custom experiment files.<BR>
