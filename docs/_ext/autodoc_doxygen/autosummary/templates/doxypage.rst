.. _{{ name }}:
{# comment
When the name is provided, we get "(INFO/1) Duplicate implicit target name:"
without this, we get undefined reference.  This needs to be fixed later.
#}

{{ underline }}
{{ title }}
{{ underline }}

{% for line in text %}
{{ line }}
{% endfor %}
{% if footnotes %}

.. rubric:: Footnotes

{% for line in footnotes %}
.. [#] {{ line }}
{% endfor %}
{% endif %}
