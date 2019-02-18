c = get_config()

c.Exporter.template_file = "./templates/nb.tpl"
c.TagRemovePreprocessor.remove_input_tags.add("rm_in")
c.TagRemovePreprocessor.remove_all_outputs_tags.add("rm_out")
