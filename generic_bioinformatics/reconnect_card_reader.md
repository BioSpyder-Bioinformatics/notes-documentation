# Reconnect card reader

The license for the card reader lives in a VM in the server clin1
- ip 192.168.168.40
- username gioele
- password biospyder

List all the virtual machines in the server
- (sudo) virsh list --all
The machine you want to restart is 'knime'
- (sudo) virsh start knime
Check if knime is up by pinging it
- ping 192.168.168.48

Now you need to restart the paxton server with the VM IP


