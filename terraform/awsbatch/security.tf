//security.tf

resource "aws_security_group" "allow_ssh" {

  name        = "allow_ssh-${random_string.rand.result}"
  description = "default ssh (22) access with Terraform"
  ingress {
    cidr_blocks = [
      "0.0.0.0/0"
    ]
    from_port = 22
    to_port   = 22
    protocol  = "tcp"
  }
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

resource "aws_security_group" "allow_http" {

  name        = "allow_http-${random_string.rand.result}"
  description = "default HTTP (80) access with Terraform"
  ingress {
    cidr_blocks = [
      "0.0.0.0/0"
    ]
    from_port = 80
    to_port   = 80
    protocol  = "tcp"
  }
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

resource "aws_security_group" "allow_shiny" {

  name        = "allow_shiny-${random_string.rand.result}"
  description = "default shiny (3838) access with Terraform"
  ingress {
    cidr_blocks = [
      "0.0.0.0/0"
    ]
    from_port = 3838
    to_port   = 3838
    protocol  = "tcp"
  }
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

resource "aws_security_group" "allow_all" {

  name        = "allow_all-${random_string.rand.result}"
  description = "default VPC security group with Terraform"
  ingress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

// Role for the cluster
resource "aws_iam_role" "ClusterRole" {
  name = "ClusterRole-${random_string.rand.result}"
  assume_role_policy = jsonencode({
    "Version" : "2012-10-17",
    "Statement" : [
      {
        "Effect" : "Allow",
        "Principal" : {
          "Service" : "batch.amazonaws.com"
        },
        "Action" : "sts:AssumeRole"
      },
      {
        "Effect" : "Allow",
        "Principal" : {
          "Service" : "ec2.amazonaws.com"
        },
        "Action" : "sts:AssumeRole"
      }
    ]
  })

}

// Role for the cluster nodes
resource "aws_iam_role" "ComputeInstanceRole" {
  name = "ComputeInstanceRole-${random_string.rand.result}"
  assume_role_policy = jsonencode({
    "Version" : "2012-10-17",
    "Statement" : [
      {
        "Effect" : "Allow",
        "Principal" : {
          "Service" : "ec2.amazonaws.com"
        },
        "Action" : "sts:AssumeRole"
      }
    ]
  })

}

// Additional role for fleeting cluster nodes
resource "aws_iam_role" "ClusterFleetRole" {
  name = "ClusterFleetRole-${random_string.rand.result}"
  assume_role_policy = jsonencode({
    "Version" : "2012-10-17",
    "Statement" : [
      {
        "Sid" : "",
        "Effect" : "Allow",
        "Principal" : {
          "Service" : "spotfleet.amazonaws.com"
        },
        "Action" : "sts:AssumeRole"
      },
      {
        "Effect" : "Allow",
        "Principal" : {
          "Service" : "ec2.amazonaws.com"
        },
        "Action" : "sts:AssumeRole"
      }
    ]
  })

}

// Role for the entrypoint
resource "aws_iam_role" "Multiaccess" {
  name = "Multiaccess-${random_string.rand.result}"
  assume_role_policy = jsonencode({
    "Version" : "2012-10-17",
    "Statement" : [
      {
        "Effect" : "Allow",
        "Principal" : {
          "Service" : "batch.amazonaws.com"
        },
        "Action" : "sts:AssumeRole"
      },
      {
        "Effect" : "Allow",
        "Principal" : {
          "Service" : "ec2.amazonaws.com"
        },
        "Action" : "sts:AssumeRole"
      }
    ]
  })

}

resource "aws_iam_instance_profile" "ComputeInstanceProfile" {
  name = "ComputeInstanceProfile-${random_string.rand.result}"
  role = aws_iam_role.ComputeInstanceRole.name
}

resource "aws_iam_instance_profile" "Multiprofile" {
  name = "Multiprofile-${random_string.rand.result}"
  role = aws_iam_role.Multiaccess.name
}

resource "aws_iam_policy_attachment" "AWSBatchServiceRole-policy-attachment" {

  name       = "AWSBatchServiceRole-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole"
  groups     = []
  users      = []
  roles      = [aws_iam_role.ClusterRole.name, aws_iam_role.Multiaccess.name]

}

resource "aws_iam_policy_attachment" "AWSTransferLoggingAccess-policy-attachment" {

  name       = "AWSTransferLoggingAccess-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSTransferLoggingAccess"
  groups     = []
  users      = []
  roles      = [aws_iam_role.ClusterRole.name]

}

resource "aws_iam_policy_attachment" "AmazonEC2ContainerServiceAutoscaleRole-policy-attachment" {

  name       = "AmazonEC2ContainerServiceAutoscaleRole-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceAutoscaleRole"
  groups     = []
  users      = []
  roles      = [aws_iam_role.ClusterRole.name]

}

resource "aws_iam_policy_attachment" "CloudWatchLogsFullAccess-policy-attachment" {

  name       = "CloudWatchLogsFullAccess-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/CloudWatchLogsFullAccess"
  groups     = []
  users      = []
  roles      = [aws_iam_role.ClusterRole.name, aws_iam_role.Multiaccess.name]

}

resource "aws_iam_policy_attachment" "AmazonEC2ContainerServiceforEC2Role-policy-attachment" {

  name       = "AmazonEC2ContainerServiceforEC2Role-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
  groups     = []
  users      = []
  roles      = [aws_iam_role.ClusterRole.name, aws_iam_role.ComputeInstanceRole.name]

}

resource "aws_iam_policy_attachment" "AmazonEC2ContainerServiceRole-policy-attachment" {

  name       = "AmazonEC2ContainerServiceRole-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceRole"
  groups     = []
  users      = []
  roles      = [aws_iam_role.ClusterRole.name]

}

resource "aws_iam_policy_attachment" "AmazonEC2ContainerRegistryFullAccess-policy-attachment" {

  name       = "AmazonEC2ContainerRegistryFullAccess-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryFullAccess"
  groups     = []
  users      = []
  roles      = [aws_iam_role.ClusterRole.name]

}

resource "aws_iam_policy_attachment" "AmazonECS_FullAccess-policy-attachment" {

  name       = "AmazonECS_FullAccess-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/AmazonECS_FullAccess"
  groups     = []
  users      = []
  roles      = [aws_iam_role.ClusterRole.name]

}


resource "aws_iam_policy_attachment" "AmazonEC2FullAccess-policy-attachment" {

  name       = "AmazonEC2FullAccess-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/AmazonEC2FullAccess"
  groups     = []
  users      = []
  roles      = [aws_iam_role.Multiaccess.name]

}

resource "aws_iam_policy_attachment" "AWSBatchServiceEventTargetRole-policy-attachment" {

  name       = "AWSBatchServiceEventTargetRole-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSBatchServiceEventTargetRole"
  groups     = []
  users      = []
  roles      = [aws_iam_role.Multiaccess.name]

}

resource "aws_iam_policy_attachment" "AmazonS3FullAccess-policy-attachment" {

  name       = "AmazonS3FullAccess-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/AmazonS3FullAccess"
  groups     = []
  users      = []
  roles      = [aws_iam_role.ClusterRole.name, aws_iam_role.Multiaccess.name, aws_iam_role.ComputeInstanceRole.name]

}

resource "aws_iam_policy_attachment" "CloudWatchFullAccess-policy-attachment" {

  name       = "CloudWatchFullAccess-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/CloudWatchFullAccess"
  groups     = []
  users      = []
  roles      = [aws_iam_role.Multiaccess.name]

}

resource "aws_iam_policy_attachment" "AmazonSSMManagedInstanceCore-policy-attachment" {

  name       = "AmazonSSMManagedInstanceCore-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore"
  groups     = []
  users      = []
  roles      = [aws_iam_role.Multiaccess.name]

}

resource "aws_iam_policy_attachment" "AWSBatchFullAccess-policy-attachment" {

  name       = "AWSBatchFullAccess-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/AWSBatchFullAccess"
  groups     = []
  users      = []
  roles      = [aws_iam_role.Multiaccess.name]

}

resource "aws_iam_policy_attachment" "CloudWatchEventsFullAccess-policy-attachment" {

  name       = "CloudWatchEventsFullAccess-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/CloudWatchEventsFullAccess"
  groups     = []
  users      = []
  roles      = [aws_iam_role.Multiaccess.name]

}

resource "aws_iam_policy_attachment" "AmazonEC2SpotFleetTaggingRole-policy-attachment" {

  name       = "AmazonEC2SpotFleetTaggingRole-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetTaggingRole"
  groups     = []
  users      = []
  roles      = [aws_iam_role.ClusterFleetRole.name]

}

resource "aws_iam_policy_attachment" "AmazonEC2SpotFleetAutoscaleRole-policy-attachment" {

  name       = "AmazonEC2SpotFleetAutoscaleRole-policy-attachment-${random_string.rand.result}"
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetAutoscaleRole"
  groups     = []
  users      = []
  roles      = [aws_iam_role.ClusterFleetRole.name]

}



// THIS MAY BE NEEDED. Spot service_linked_role
// resource "aws_iam_service_linked_role" "spot_service_linked_role" {
//   aws_service_name = "spot.amazonaws.com"
// }
