#!/bin/bash

# Let's update first
sudo yum update -y

sudo mkdir -p /mnt/${bucket_prefix}-${count}

sudo s3fs -o iam_role="Multiaccess-${rand}" -o url="https://s3-${region}.amazonaws.com/" -o endpoint=${region} -o dbglevel=info -o umask=0022 -o uid=1000 -o gid=1000 -o curldbg -o allow_other -o default_acl=${bucket_acl} -o use_cache=/tmp ${bucket_prefix}-${count} /mnt/${bucket_prefix}-${count}

sudo sed -i '/^PasswordAuthentication/c\PasswordAuthentication yes' /etc/ssh/sshd_config

sudo echo "ec2-user:${ec2_password}"|chpasswd

sudo systemctl restart sshd

mkdir -p /home/ec2-user/git

# Git of MoP2

cd /home/ec2-user/git; git clone --recurse-submodules ${repourl}

sudo chown -R ec2-user:ec2-user /home/ec2-user/git

# This is for Singularity

sudo yum install -y debootstrap

# Let's record .bash_history for activity tracking

cat <<EOF >> /home/ec2-user/.bashrc
HISTFILESIZE=400000000
HISTSIZE=10000
PROMPT_COMMAND="history -a"
shopt -s histappend
EOF

# We clean above history
rm /home/ec2-user/.bash_history; touch /home/ec2-user/.bash_history; chown ec2-user:ec2-user /home/ec2-user/.bash_history;
