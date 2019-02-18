c = get_config()

c.Exporter.template_file = "./templates/slides.tpl"
c.TagRemovePreprocessor.remove_input_tags.add("rm_in")
c.TagRemovePreprocessor.remove_all_outputs_tags.add("rm_out")
c.SlidesExporter.reveal_theme = "simple"
exporter_settings = {"exclude_input_prompt": True, "exclude_output_prompt": True}
c.TemplateExporter.update(exporter_settings)
