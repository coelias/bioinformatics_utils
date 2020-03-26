#!/usr/bin/python
import subprocess
import re
from os.path import expanduser

default='''IdentityFile ~/.ssh/id_rsa
ForwardAgent yes



'''

#n1=0
#a=subprocess.Popen(['ssh','bastion-dev','./list_servers.sh'],stdout=subprocess.PIPE)
#for i in a.stdout.read().split('\n'):
#    n1+=1
#    i=i.strip()
#    if i:
#        _tmp=re.findall('\S+',i)
#        name="_".join(_tmp[:-1])
#        ip=_tmp[-1]
#        default+="Host PROD_{0}\n\tHostName {1}\n\tUser ubuntu\n\tProxyCommand ssh ubuntu@bastion-dev nc %h %p\n\n".format(name,ip)
#
#n2=0
#a=subprocess.Popen(['ssh','bastion-prod','./list_servers.sh'],stdout=subprocess.PIPE)
#for i in a.stdout.read().split('\n'):
#    n2+=1
#    i=i.strip()
#    if i:
#        _tmp=re.findall('\S+',i)
#        name="_".join(_tmp[:-1])
#        ip=_tmp[-1]
#        default+="Host PROD_{0}\n\tHostName {1}\n\tUser ubuntu\n\tProxyCommand  ssh ubuntu@bastion-prod nc %h %p\n\n".format(name,ip)
#
#a=subprocess.Popen(['ssh','self-service','./list_instances.sh'],stdout=subprocess.PIPE)
#a=[i.strip() for i in a.stdout if 'uid:' in i or i.startswith('10.')]
#n3=len(a)
#ipname=zip(a[::2],a[1::2])
#
#for ip,name in ipname:
#    name=name[5:]
#    default+="Host SELF_{0}\n\tHostName {1}\n\tUser ubuntu\n\n".format(name,ip)
#
#
#if all([n1,n2,n3]):
#    with open(expanduser('~/.ssh/config'),'w') as f:
#        f.write(default)
