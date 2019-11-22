//batch.tf

resource "aws_batch_compute_environment" "nf-spot" {
  
    compute_environment_name = "nf-spot"

    compute_resources {
      // instance_role = "${aws_iam_instance_profile.S3access.arn}"
      instance_role = "arn:aws:iam::132458770246:instance-profile/S3access"
      bid_percentage = 50
  
      image_id = "ami-0005d00d062a1f312"
      //launch_template  {
      //  launch_template_id = "${aws_launch_template.test.id}"
      //  version = "$Latest"
      //}
  
      max_vcpus = 16
      min_vcpus = 0
      desired_vcpus = 0
  
      instance_type = [ "optimal" ]
  
      subnets = ["subnet-8a280df7", "subnet-c54d6588", "subnet-b85ab5d2"]
  
      spot_iam_fleet_role = "arn:aws:iam::132458770246:role/AmazonEC2SpotFleetRole"
  
      type = "SPOT"
          
      security_group_ids = [
        "sg-c290cbac"
      ]
      
    }

    // service_role = "${aws_iam_role.AWSBatchServiceRole.arn}"
    service_role = "arn:aws:iam::132458770246:role/service-role/AWSBatchServiceRole"
    type         = "MANAGED"
    depends_on   = ["aws_iam_policy_attachment.AWSBatchServiceRole-policy-attachment"]
    

}

resource "aws_batch_job_queue" "spot" {
  name                 = "spot"
  state                = "ENABLED"
  priority             = 1
  compute_environments = ["${aws_batch_compute_environment.nf-spot.arn}"]
  
}


